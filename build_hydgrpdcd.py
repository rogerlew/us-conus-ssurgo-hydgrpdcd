#!/usr/bin/env python3
"""Generate a CONUS hydrologic soil group raster from gNATSGO mukey data.

The output raster encodes canonical HSG groups using:
  - 0: NoData / unresolved
  - 1: A
  - 2: B
  - 3: C
  - 4: D

Data flow:
1) Read unique mukeys from the gNATSGO VAT DBF.
2) Query NRCS SDA muaggatt.hydgrpdcd for missing mukeys (batched).
3) Optionally fall back to dominant component hydgrp when muaggatt is unresolved.
4) Canonicalize hydgrpdcd/hydgrp values to A|B|C|D.
5) Remap source mukey raster block-wise into coded HSG raster.
6) Emit lookup CSV + metadata JSON for auditability.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import requests
from osgeo import gdal, ogr
from osgeo.gdalconst import GDT_Byte, GA_ReadOnly

gdal.UseExceptions()

SDA_ENDPOINT = "https://SDMDataAccess.nrcs.usda.gov/Tabular/post.rest"
CODE_TO_HSG = {0: "NODATA", 1: "A", 2: "B", 3: "C", 4: "D"}
HSG_TO_CODE = {"A": 1, "B": 2, "C": 3, "D": 4}
SIMPLE_TEXTURE_TO_HSG = {
    "sand loam": "B",
    "loam": "B",
    "silt loam": "C",
    "clay loam": "D",
}


def _chunked(values: Sequence[int], chunk_size: int) -> Iterator[Sequence[int]]:
    if chunk_size <= 0:
        raise ValueError(f"chunk_size must be > 0, got {chunk_size}")
    for idx in range(0, len(values), chunk_size):
        yield values[idx : idx + chunk_size]


def read_vat_mukeys(vat_dbf: Path, value_field: str = "Value") -> List[int]:
    """Read unique positive mukey values from a VAT DBF file."""
    dataset = ogr.Open(str(vat_dbf), 0)
    if dataset is None:
        raise FileNotFoundError(f"Unable to open VAT DBF: {vat_dbf}")

    layer = dataset.GetLayer(0)
    layer_defn = layer.GetLayerDefn()
    field_index = layer_defn.GetFieldIndex(value_field)
    if field_index < 0:
        raise KeyError(f"Field '{value_field}' not found in {vat_dbf}")

    mukeys: List[int] = []
    for feature in layer:
        value = int(feature.GetFieldAsInteger64(field_index))
        if value > 0:
            mukeys.append(value)

    dataset = None

    if not mukeys:
        raise ValueError(f"No positive mukeys found in {vat_dbf}")

    return sorted(set(mukeys))


def load_lookup_csv(
    lookup_csv: Path,
) -> Tuple[Dict[int, Optional[str]], Dict[int, str], Dict[int, str]]:
    """Load existing mukey->hydgrpdcd_raw mappings from lookup CSV."""
    lookup: Dict[int, Optional[str]] = {}
    source_lookup: Dict[int, str] = {}
    detail_lookup: Dict[int, str] = {}
    if not lookup_csv.exists():
        return lookup, source_lookup, detail_lookup

    with lookup_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"mukey", "hydgrpdcd_raw"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{lookup_csv} missing required columns: {sorted(missing)}"
            )
        for row in reader:
            mukey = int(row["mukey"])
            raw = (row.get("hydgrpdcd_raw") or "").strip()
            lookup[mukey] = raw if raw else None
            raw_source = (row.get("raw_source") or "lookup_cache").strip()
            source_lookup[mukey] = raw_source or "lookup_cache"
            raw_detail = (row.get("raw_detail") or "").strip()
            detail_lookup[mukey] = raw_detail

    return lookup, source_lookup, detail_lookup


def query_hydgrpdcd_from_sda(
    mukeys: Sequence[int],
    *,
    batch_size: int,
    timeout_seconds: int,
    max_retries: int,
) -> Dict[int, Optional[str]]:
    """Fetch muaggatt.hydgrpdcd from SDA for the supplied mukeys."""
    if not mukeys:
        return {}

    total_batches = math.ceil(len(mukeys) / batch_size)
    session = requests.Session()
    results: Dict[int, Optional[str]] = {}

    for batch_index, chunk in enumerate(_chunked(mukeys, batch_size), start=1):
        query = (
            "SELECT mukey, hydgrpdcd "
            "FROM muaggatt "
            f"WHERE mukey IN ({','.join(str(v) for v in chunk)})"
        )
        payload = {"query": query, "format": "json"}

        last_exc: Optional[Exception] = None
        for attempt in range(1, max_retries + 1):
            try:
                response = session.post(
                    SDA_ENDPOINT, data=payload, timeout=timeout_seconds
                )
                response.raise_for_status()
                body = response.json()
                rows = body.get("Table") or []
                for row in rows:
                    mukey = int(row[0])
                    raw = row[1] if len(row) > 1 else None
                    raw_clean = str(raw).strip() if raw is not None else ""
                    results[mukey] = raw_clean if raw_clean else None
                break
            except (requests.RequestException, ValueError) as exc:
                last_exc = exc
                if attempt >= max_retries:
                    raise RuntimeError(
                        "Failed to fetch hydgrpdcd from SDA for batch "
                        f"{batch_index}/{total_batches}"
                    ) from exc
                sleep_seconds = 2 ** (attempt - 1)
                print(
                    f"[warn] SDA batch {batch_index}/{total_batches} attempt "
                    f"{attempt}/{max_retries} failed ({exc}); retrying in "
                    f"{sleep_seconds}s..."
                )
                time.sleep(sleep_seconds)
        else:
            if last_exc is not None:
                raise RuntimeError("Unreachable SDA retry state") from last_exc

        if batch_index % 10 == 0 or batch_index == total_batches:
            print(f"[info] SDA batch {batch_index}/{total_batches} complete")

    return results


def _select_component_hydgrp_rows(
    rows: Sequence[Sequence[object]],
) -> Dict[int, Optional[str]]:
    """Select dominant non-empty component hydgrp for each mukey.

    Expects rows pre-sorted by ``mukey, comppct_r DESC, cokey``.
    """
    selected: Dict[int, Optional[str]] = {}

    for row in rows:
        if len(row) < 4:
            continue
        mukey = int(row[0])
        if mukey in selected and selected[mukey] is not None:
            continue

        hydgrp = row[3]
        value = str(hydgrp).strip().upper() if hydgrp is not None else ""
        if value:
            selected[mukey] = value
        elif mukey not in selected:
            selected[mukey] = None

    return selected


def query_component_hydgrp_major_from_sda(
    mukeys: Sequence[int],
    *,
    batch_size: int,
    timeout_seconds: int,
    max_retries: int,
) -> Dict[int, Optional[str]]:
    """Fetch dominant component hydgrp values for supplied mukeys."""
    if not mukeys:
        return {}

    total_batches = math.ceil(len(mukeys) / batch_size)
    session = requests.Session()
    results: Dict[int, Optional[str]] = {}

    for batch_index, chunk in enumerate(_chunked(mukeys, batch_size), start=1):
        query = (
            "SELECT mukey, cokey, comppct_r, hydgrp "
            "FROM component "
            f"WHERE mukey IN ({','.join(str(v) for v in chunk)}) "
            "ORDER BY mukey, comppct_r DESC, cokey"
        )
        payload = {"query": query, "format": "json"}

        last_exc: Optional[Exception] = None
        for attempt in range(1, max_retries + 1):
            try:
                response = session.post(
                    SDA_ENDPOINT, data=payload, timeout=timeout_seconds
                )
                response.raise_for_status()
                body = response.json()
                rows = body.get("Table") or []
                selected = _select_component_hydgrp_rows(rows)
                results.update(selected)
                break
            except (requests.RequestException, ValueError) as exc:
                last_exc = exc
                if attempt >= max_retries:
                    raise RuntimeError(
                        "Failed to fetch component hydgrp from SDA for batch "
                        f"{batch_index}/{total_batches}"
                    ) from exc
                sleep_seconds = 2 ** (attempt - 1)
                print(
                    f"[warn] component fallback batch {batch_index}/{total_batches} "
                    f"attempt {attempt}/{max_retries} failed ({exc}); retrying "
                    f"in {sleep_seconds}s..."
                )
                time.sleep(sleep_seconds)
        else:
            if last_exc is not None:
                raise RuntimeError("Unreachable component fallback retry state") from last_exc

        if batch_index % 10 == 0 or batch_index == total_batches:
            print(
                f"[info] component fallback batch {batch_index}/{total_batches} complete"
            )

    return results


def _soil_texture_detailed(clay: float, sand: float) -> Optional[str]:
    """Return the detailed USDA soil texture class.

    Mirrors WEPPpy's ``_soil_texture`` implementation in
    ``wepppy/wepp/soils/utils/utils.py``.
    """
    assert sand + clay <= 100
    silt = 100.0 - sand - clay

    if clay >= 40:
        if silt >= 40:
            return "silty clay"
        if sand <= 45:
            return "clay"

    if clay >= 35 and sand > 45:
        return "sandy clay"

    if clay >= 27:
        if sand <= 20:
            return "silty clay loam"
        if sand <= 45:
            return "clay loam"
    else:
        if silt >= 50:
            if clay < 12.0 and silt >= 80:
                return "silt"
            return "silt loam"
        if silt >= 28 and clay >= 7 and sand <= 52:
            return "loam"

    if clay >= 20 and sand > 45 and silt <= 28:
        return "sandy clay loam"

    if silt + 1.5 * clay < 15:
        return "sand"
    if silt + 2 * clay < 30:
        return "loamy sand"
    return "sandy loam"


def simple_texture_4class(clay: float, sand: float) -> Optional[str]:
    """Return WEPPpy's 4-class simplified texture label.

    Mirrors WEPPpy's ``simple_texture`` implementation in
    ``wepppy/wepp/soils/utils/utils.py``.
    """
    cs = clay + sand
    if (clay <= 27.0 and cs <= 50.0) or (clay > 27.0 and sand <= 20.0 and cs <= 50.0):
        return "silt loam"
    if (6.0 <= clay <= 27.0) and (50.0 < cs <= 72.0) and sand <= 52:
        return "loam"
    if (sand > 52 or cs > 50 and clay < 6) and sand >= 50:
        return "sand loam"
    if (cs > 72 and sand < 50) or (clay > 27 and (20 < sand <= 45)) or (sand <= 20 and cs > 50):
        return "clay loam"

    detailed = _soil_texture_detailed(clay, sand)
    if detailed is None:
        return None
    if detailed.startswith("sand"):
        return "sand loam"
    if detailed.startswith("silt"):
        return "silt loam"
    if detailed.startswith("clay"):
        return "clay loam"
    if detailed.startswith("loam"):
        return "loam"
    return None


def _to_float_or_none(value: object) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def map_horizon_simple_texture_to_hsg(
    *,
    sand_value: object,
    clay_value: object,
) -> Tuple[Optional[str], str, str]:
    """Map horizon sand/clay to HSG using WEPPpy 4-class simple texture."""
    sand = _to_float_or_none(sand_value)
    clay = _to_float_or_none(clay_value)
    if sand is None or clay is None:
        return None, "texture4_missing_sand_or_clay", ""

    if sand < 0 or clay < 0 or sand > 100 or clay > 100 or (sand + clay) > 100:
        return None, "texture4_invalid_sand_clay", ""

    texture4 = simple_texture_4class(clay=clay, sand=sand)
    if texture4 is None:
        return None, "texture4_unclassified", ""

    hsg = SIMPLE_TEXTURE_TO_HSG.get(texture4)
    if hsg is None:
        return None, "texture4_unmapped_class", texture4

    return hsg, "texture4_mapped", texture4


def _select_chorizon_texture_fallback_rows(
    rows: Sequence[Sequence[object]],
) -> Dict[int, Tuple[Optional[str], str, str]]:
    """Select per-mukey HSG from ordered component/horizon rows."""
    selected: Dict[int, Tuple[Optional[str], str, str]] = {}

    for row in rows:
        if len(row) < 9:
            continue
        mukey = int(row[0])
        if mukey in selected and selected[mukey][0] is not None:
            continue

        hsg, status, texture4 = map_horizon_simple_texture_to_hsg(
            sand_value=row[6],
            clay_value=row[7],
        )
        raw_texture_text = str(row[8]).strip() if row[8] is not None else ""
        detail = texture4
        if raw_texture_text:
            detail = f"{texture4}|chtexturegrp={raw_texture_text}" if texture4 else raw_texture_text
        if hsg is not None:
            selected[mukey] = (hsg, status, detail)
        else:
            selected[mukey] = (None, status, detail)

    return selected


def query_chorizon_chtexturegrp_fallback_from_sda(
    mukeys: Sequence[int],
    *,
    batch_size: int,
    timeout_seconds: int,
    max_retries: int,
) -> Tuple[Dict[int, Optional[str]], Dict[int, str], Counter]:
    """Fetch HSG fallback values from ordered component/chorizon/chtexture rows."""
    if not mukeys:
        return {}, {}, Counter()

    total_batches = math.ceil(len(mukeys) / batch_size)
    session = requests.Session()
    hsg_lookup: Dict[int, Optional[str]] = {}
    detail_lookup: Dict[int, str] = {}
    status_counts: Counter = Counter()

    for batch_index, chunk in enumerate(_chunked(mukeys, batch_size), start=1):
        query = (
            "SELECT c.mukey, c.cokey, c.comppct_r, h.chkey, h.hzdept_r, h.hzdepb_r, "
            "h.sandtotal_r, h.claytotal_r, t.texture "
            "FROM component c "
            "INNER JOIN chorizon h ON c.cokey = h.cokey "
            "LEFT JOIN chtexturegrp t ON h.chkey = t.chkey "
            f"WHERE c.mukey IN ({','.join(str(v) for v in chunk)}) "
            "ORDER BY c.mukey, c.comppct_r DESC, h.hzdept_r ASC, h.hzdepb_r ASC, h.chkey"
        )
        payload = {"query": query, "format": "json"}

        last_exc: Optional[Exception] = None
        for attempt in range(1, max_retries + 1):
            try:
                response = session.post(
                    SDA_ENDPOINT, data=payload, timeout=timeout_seconds
                )
                response.raise_for_status()
                body = response.json()
                rows = body.get("Table") or []
                selected = _select_chorizon_texture_fallback_rows(rows)
                for mukey, (hsg, status, raw_texture) in selected.items():
                    hsg_lookup[mukey] = hsg
                    detail_lookup[mukey] = raw_texture
                    status_counts[status] += 1
                break
            except (requests.RequestException, ValueError) as exc:
                last_exc = exc
                if attempt >= max_retries:
                    raise RuntimeError(
                        "Failed to fetch chorizon/chtexturegrp fallback from SDA "
                        f"for batch {batch_index}/{total_batches}"
                    ) from exc
                sleep_seconds = 2 ** (attempt - 1)
                print(
                    f"[warn] chorizon fallback batch {batch_index}/{total_batches} "
                    f"attempt {attempt}/{max_retries} failed ({exc}); retrying "
                    f"in {sleep_seconds}s..."
                )
                time.sleep(sleep_seconds)
        else:
            if last_exc is not None:
                raise RuntimeError("Unreachable chorizon fallback retry state") from last_exc

        for mukey in chunk:
            if mukey not in hsg_lookup:
                hsg_lookup[mukey] = None
                detail_lookup[mukey] = ""
                status_counts["texture_no_rows"] += 1

        if batch_index % 10 == 0 or batch_index == total_batches:
            print(
                f"[info] chorizon fallback batch {batch_index}/{total_batches} complete"
            )

    return hsg_lookup, detail_lookup, status_counts


def canonicalize_hydgrpdcd(
    raw_value: Optional[str], *, dual_group_policy: str
) -> Tuple[Optional[str], str]:
    """Normalize raw hydgrpdcd values into canonical A|B|C|D or unresolved."""
    if raw_value is None:
        return None, "missing_raw"

    value = str(raw_value).strip().upper()
    if not value:
        return None, "missing_raw"

    if value in HSG_TO_CODE:
        return value, "ok"

    if value in {"A/D", "B/D", "C/D"}:
        if dual_group_policy == "assume_d":
            return "D", "ok_dual_assume_d"
        if dual_group_policy == "error":
            return None, "dual_group_blocked"
        raise ValueError(f"Unknown dual_group_policy: {dual_group_policy}")

    return None, "unknown_raw_code"


def build_lookup_records(
    mukeys: Sequence[int],
    raw_lookup: Dict[int, Optional[str]],
    *,
    dual_group_policy: str,
    raw_source_lookup: Optional[Dict[int, str]] = None,
    raw_detail_lookup: Optional[Dict[int, str]] = None,
) -> Tuple[np.ndarray, List[Dict[str, object]], Counter, Counter]:
    """Build mukey LUT and audit records from raw hydgrpdcd lookups."""
    if not mukeys:
        raise ValueError("mukeys cannot be empty")

    lut = np.zeros(max(mukeys) + 1, dtype=np.uint8)
    records: List[Dict[str, object]] = []
    status_counts: Counter = Counter()
    source_counts: Counter = Counter()

    for mukey in mukeys:
        raw = raw_lookup.get(mukey)
        raw_source = (
            raw_source_lookup.get(mukey, "unknown")
            if raw_source_lookup is not None
            else "unknown"
        )
        raw_detail = (
            raw_detail_lookup.get(mukey, "")
            if raw_detail_lookup is not None
            else ""
        )
        canonical, status = canonicalize_hydgrpdcd(
            raw, dual_group_policy=dual_group_policy
        )
        code = HSG_TO_CODE.get(canonical, 0)
        lut[mukey] = code

        records.append(
            {
                "mukey": mukey,
                "hydgrpdcd_raw": raw or "",
                "hsg_group": canonical or "",
                "code": code,
                "status": status,
                "raw_source": raw_source,
                "raw_detail": raw_detail,
            }
        )
        status_counts[status] += 1
        source_counts[raw_source] += 1

    return lut, records, status_counts, source_counts


def write_lookup_csv(records: Iterable[Dict[str, object]], lookup_csv: Path) -> None:
    lookup_csv.parent.mkdir(parents=True, exist_ok=True)
    with lookup_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mukey",
                "hydgrpdcd_raw",
                "hsg_group",
                "code",
                "status",
                "raw_source",
                "raw_detail",
            ],
        )
        writer.writeheader()
        for row in records:
            writer.writerow(row)


def write_coded_raster(
    source_raster: Path,
    output_raster: Path,
    lut: np.ndarray,
    *,
    nodata_code: int = 0,
) -> Dict[str, object]:
    """Write output coded raster by remapping source mukey values using LUT."""
    src_ds = gdal.Open(str(source_raster), GA_ReadOnly)
    if src_ds is None:
        raise FileNotFoundError(f"Unable to open source raster: {source_raster}")

    src_band = src_ds.GetRasterBand(1)
    x_size = src_ds.RasterXSize
    y_size = src_ds.RasterYSize
    block_x, block_y = src_band.GetBlockSize()
    if block_x <= 0 or block_y <= 0:
        block_x, block_y = 1024, 1024

    output_raster.parent.mkdir(parents=True, exist_ok=True)

    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(
        str(output_raster),
        x_size,
        y_size,
        1,
        GDT_Byte,
        options=[
            "TILED=YES",
            "COMPRESS=LZW",
            "BIGTIFF=IF_SAFER",
        ],
    )
    if dst_ds is None:
        raise RuntimeError(f"Failed to create output raster: {output_raster}")

    dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
    dst_ds.SetProjection(src_ds.GetProjection())
    dst_band = dst_ds.GetRasterBand(1)
    dst_band.SetNoDataValue(nodata_code)

    cell_counts = np.zeros(max(CODE_TO_HSG.keys()) + 1, dtype=np.int64)
    lookup_max = lut.shape[0]
    total_blocks = math.ceil(y_size / block_y)

    for block_idx, y_off in enumerate(range(0, y_size, block_y), start=1):
        rows = min(block_y, y_size - y_off)
        src_arr = src_band.ReadAsArray(0, y_off, x_size, rows)
        if src_arr is None:
            raise RuntimeError(f"Failed reading source block at y={y_off}")

        out_arr = np.zeros(src_arr.shape, dtype=np.uint8)
        valid = (src_arr >= 0) & (src_arr < lookup_max)
        out_arr[valid] = lut[src_arr[valid]]

        dst_band.WriteArray(out_arr, 0, y_off)
        cell_counts += np.bincount(
            out_arr.ravel(), minlength=len(cell_counts)
        ).astype(np.int64)

        if block_idx % 100 == 0 or block_idx == total_blocks:
            print(f"[info] raster block row {block_idx}/{total_blocks} complete")

    dst_band.FlushCache()
    dst_ds.FlushCache()
    dst_ds = None
    src_ds = None

    return {
        "width": x_size,
        "height": y_size,
        "cell_counts_by_code": {
            str(code): int(cell_counts[code]) for code in range(len(cell_counts))
        },
    }


def write_metadata_json(
    metadata_path: Path,
    *,
    source_raster: Path,
    source_vat: Path,
    output_raster: Path,
    lookup_csv: Path,
    dual_group_policy: str,
    status_counts: Counter,
    source_counts: Counter,
    raster_summary: Dict[str, object],
    mukey_count: int,
    component_fallback_enabled: bool,
    unresolved_before_component_fallback: int,
    unresolved_after_component_fallback: int,
    component_fallback_applied: int,
    chorizon_fallback_enabled: bool,
    unresolved_before_chorizon_fallback: int,
    unresolved_after_chorizon_fallback: int,
    chorizon_fallback_applied: int,
    chorizon_texture_status_counts: Counter,
) -> None:
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_raster": str(source_raster),
        "source_vat": str(source_vat),
        "output_raster": str(output_raster),
        "lookup_csv": str(lookup_csv),
        "sda_endpoint": SDA_ENDPOINT,
        "dual_group_policy": dual_group_policy,
        "codebook": CODE_TO_HSG,
        "mukey_count": mukey_count,
        "lookup_status_counts": {
            key: int(value) for key, value in sorted(status_counts.items())
        },
        "lookup_raw_source_counts": {
            key: int(value) for key, value in sorted(source_counts.items())
        },
        "component_fallback_enabled": component_fallback_enabled,
        "unresolved_before_component_fallback": unresolved_before_component_fallback,
        "unresolved_after_component_fallback": unresolved_after_component_fallback,
        "component_fallback_applied": component_fallback_applied,
        "chorizon_fallback_enabled": chorizon_fallback_enabled,
        "unresolved_before_chorizon_fallback": unresolved_before_chorizon_fallback,
        "unresolved_after_chorizon_fallback": unresolved_after_chorizon_fallback,
        "chorizon_fallback_applied": chorizon_fallback_applied,
        "chorizon_texture_status_counts": {
            key: int(value)
            for key, value in sorted(chorizon_texture_status_counts.items())
        },
        "texture4_to_hsg_fallback_map": SIMPLE_TEXTURE_TO_HSG,
        "raster_summary": raster_summary,
    }
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build hydgrpdcd coded raster from gNATSGO mukey raster."
    )
    parser.add_argument(
        "--source-raster",
        type=Path,
        default=Path("/home/geodata/ssurgo/gNATSGSO/2025/gNATSGO_mukey_202502.tif"),
        help="Input mukey GeoTIFF.",
    )
    parser.add_argument(
        "--source-vat",
        type=Path,
        default=Path(
            "/home/geodata/ssurgo/gNATSGSO/2025/gNATSGO_mukey_202502.tif.vat.dbf"
        ),
        help="Input VAT DBF containing mukey values.",
    )
    parser.add_argument(
        "--output-raster",
        type=Path,
        default=Path("/home/geodata/ssurgo/hydgrpdcd/hydgrpdcd.tif"),
        help="Output HSG-coded GeoTIFF.",
    )
    parser.add_argument(
        "--lookup-csv",
        type=Path,
        default=Path("/home/geodata/ssurgo/hydgrpdcd/hydgrpdcd_lookup.csv"),
        help="Lookup cache CSV (mukey -> raw hydgrpdcd -> canonical code).",
    )
    parser.add_argument(
        "--metadata-json",
        type=Path,
        default=Path("/home/geodata/ssurgo/hydgrpdcd/hydgrpdcd_metadata.json"),
        help="Metadata output JSON path.",
    )
    parser.add_argument(
        "--refresh-lookup",
        action="store_true",
        help="Ignore existing lookup CSV and refetch all mukeys from SDA.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=2000,
        help="SDA query mukey batch size.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=60,
        help="Per-request timeout for SDA calls.",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=4,
        help="Retry attempts per SDA batch request.",
    )
    parser.add_argument(
        "--dual-group-policy",
        choices=("assume_d", "error"),
        default="assume_d",
        help=(
            "How to handle dual classes (A/D, B/D, C/D): "
            "'assume_d' maps to D; 'error' leaves unresolved (NoData)."
        ),
    )
    parser.add_argument(
        "--component-fallback",
        dest="component_fallback",
        action="store_true",
        default=True,
        help=(
            "Fallback unresolved muaggatt hydgrpdcd values to dominant "
            "component hydgrp."
        ),
    )
    parser.add_argument(
        "--no-component-fallback",
        dest="component_fallback",
        action="store_false",
        help="Disable dominant-component hydgrp fallback.",
    )
    parser.add_argument(
        "--chorizon-fallback",
        dest="chorizon_fallback",
        action="store_true",
        default=True,
        help=(
            "Fallback unresolved values using ordered component/chorizon/"
            "chtexturegrp texture mapping."
        ),
    )
    parser.add_argument(
        "--no-chorizon-fallback",
        dest="chorizon_fallback",
        action="store_false",
        help="Disable chorizon/chtexturegrp fallback.",
    )
    return parser


def unresolved_mukeys(
    mukeys: Sequence[int],
    raw_lookup: Dict[int, Optional[str]],
    *,
    dual_group_policy: str,
) -> List[int]:
    """Return mukeys whose current raw values do not canonicalize."""
    unresolved: List[int] = []
    for mukey in mukeys:
        canonical, _status = canonicalize_hydgrpdcd(
            raw_lookup.get(mukey), dual_group_policy=dual_group_policy
        )
        if canonical is None:
            unresolved.append(mukey)
    return unresolved


def main() -> int:
    args = build_parser().parse_args()

    print(f"[info] reading mukeys from {args.source_vat}")
    mukeys = read_vat_mukeys(args.source_vat)
    print(f"[info] unique mukeys in VAT: {len(mukeys)}")

    raw_lookup: Dict[int, Optional[str]] = {}
    raw_source_lookup: Dict[int, str] = {}
    raw_detail_lookup: Dict[int, str] = {}
    if args.lookup_csv.exists() and not args.refresh_lookup:
        raw_lookup, raw_source_lookup, raw_detail_lookup = load_lookup_csv(
            args.lookup_csv
        )
        print(
            f"[info] loaded lookup cache rows: {len(raw_lookup)} from {args.lookup_csv}"
        )

    missing = [mukey for mukey in mukeys if mukey not in raw_lookup]
    if missing:
        print(f"[info] querying SDA for {len(missing)} missing mukeys")
        fetched = query_hydgrpdcd_from_sda(
            missing,
            batch_size=args.batch_size,
            timeout_seconds=args.timeout_seconds,
            max_retries=args.max_retries,
        )
        raw_lookup.update(fetched)
        for mukey in fetched:
            raw_source_lookup[mukey] = "muaggatt"
            raw_detail_lookup.setdefault(mukey, "")
        # Preserve absent rows as explicit missing values for auditability.
        for mukey in missing:
            raw_lookup.setdefault(mukey, None)
            raw_source_lookup.setdefault(mukey, "muaggatt")
            raw_detail_lookup.setdefault(mukey, "")
    else:
        print("[info] lookup cache already covers all mukeys")

    unresolved_before_component_fallback = len(
        unresolved_mukeys(
            mukeys,
            raw_lookup,
            dual_group_policy=args.dual_group_policy,
        )
    )
    component_fallback_applied = 0
    unresolved_after_component_fallback = unresolved_before_component_fallback

    if args.component_fallback and unresolved_before_component_fallback > 0:
        unresolved_for_fallback = unresolved_mukeys(
            mukeys,
            raw_lookup,
            dual_group_policy=args.dual_group_policy,
        )
        print(
            "[info] component fallback query for "
            f"{len(unresolved_for_fallback)} unresolved mukeys"
        )
        component_fallback = query_component_hydgrp_major_from_sda(
            unresolved_for_fallback,
            batch_size=args.batch_size,
            timeout_seconds=args.timeout_seconds,
            max_retries=args.max_retries,
        )
        for mukey, raw in component_fallback.items():
            if raw is None:
                continue
            raw_lookup[mukey] = raw
            raw_source_lookup[mukey] = "component_hydgrp_major"
            raw_detail_lookup[mukey] = ""
            component_fallback_applied += 1

        unresolved_after_component_fallback = len(
            unresolved_mukeys(
                mukeys,
                raw_lookup,
                dual_group_policy=args.dual_group_policy,
            )
        )
        print(
            "[info] component fallback applied to "
            f"{component_fallback_applied} mukeys; unresolved now "
            f"{unresolved_after_component_fallback}"
        )

    unresolved_before_chorizon_fallback = unresolved_after_component_fallback
    unresolved_after_chorizon_fallback = unresolved_before_chorizon_fallback
    chorizon_fallback_applied = 0
    chorizon_texture_status_counts: Counter = Counter()

    if args.chorizon_fallback and unresolved_before_chorizon_fallback > 0:
        unresolved_for_texture = unresolved_mukeys(
            mukeys,
            raw_lookup,
            dual_group_policy=args.dual_group_policy,
        )
        print(
            "[info] chorizon/chtexturegrp fallback query for "
            f"{len(unresolved_for_texture)} unresolved mukeys"
        )
        texture_hsg, texture_detail, texture_status_counts = (
            query_chorizon_chtexturegrp_fallback_from_sda(
                unresolved_for_texture,
                batch_size=args.batch_size,
                timeout_seconds=args.timeout_seconds,
                max_retries=args.max_retries,
            )
        )
        chorizon_texture_status_counts.update(texture_status_counts)
        for mukey, hsg in texture_hsg.items():
            if hsg is None:
                continue
            raw_lookup[mukey] = hsg
            raw_source_lookup[mukey] = "chorizon_chtexturegrp_fallback"
            raw_detail_lookup[mukey] = texture_detail.get(mukey, "")
            chorizon_fallback_applied += 1

        unresolved_after_chorizon_fallback = len(
            unresolved_mukeys(
                mukeys,
                raw_lookup,
                dual_group_policy=args.dual_group_policy,
            )
        )
        print(
            "[info] chorizon fallback applied to "
            f"{chorizon_fallback_applied} mukeys; unresolved now "
            f"{unresolved_after_chorizon_fallback}"
        )

    print("[info] canonicalizing lookup and building LUT")
    lut, records, status_counts, source_counts = build_lookup_records(
        mukeys,
        raw_lookup,
        dual_group_policy=args.dual_group_policy,
        raw_source_lookup=raw_source_lookup,
        raw_detail_lookup=raw_detail_lookup,
    )
    write_lookup_csv(records, args.lookup_csv)
    print(f"[info] wrote lookup CSV: {args.lookup_csv}")

    print(f"[info] writing output raster: {args.output_raster}")
    raster_summary = write_coded_raster(
        source_raster=args.source_raster,
        output_raster=args.output_raster,
        lut=lut,
        nodata_code=0,
    )

    write_metadata_json(
        args.metadata_json,
        source_raster=args.source_raster,
        source_vat=args.source_vat,
        output_raster=args.output_raster,
        lookup_csv=args.lookup_csv,
        dual_group_policy=args.dual_group_policy,
        status_counts=status_counts,
        source_counts=source_counts,
        raster_summary=raster_summary,
        mukey_count=len(mukeys),
        component_fallback_enabled=bool(args.component_fallback),
        unresolved_before_component_fallback=unresolved_before_component_fallback,
        unresolved_after_component_fallback=unresolved_after_component_fallback,
        component_fallback_applied=component_fallback_applied,
        chorizon_fallback_enabled=bool(args.chorizon_fallback),
        unresolved_before_chorizon_fallback=unresolved_before_chorizon_fallback,
        unresolved_after_chorizon_fallback=unresolved_after_chorizon_fallback,
        chorizon_fallback_applied=chorizon_fallback_applied,
        chorizon_texture_status_counts=chorizon_texture_status_counts,
    )
    print(f"[info] wrote metadata JSON: {args.metadata_json}")
    print("[info] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

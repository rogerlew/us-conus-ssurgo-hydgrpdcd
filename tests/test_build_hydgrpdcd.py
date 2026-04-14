from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
from osgeo import gdal, osr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_hydgrpdcd as mod


def _create_test_uint32_raster(path: Path, array: np.ndarray) -> None:
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(
        str(path),
        array.shape[1],
        array.shape[0],
        1,
        gdal.GDT_UInt32,
        options=["TILED=YES", "COMPRESS=LZW"],
    )
    ds.SetGeoTransform((100.0, 30.0, 0.0, 200.0, 0.0, -30.0))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(5070)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.WriteArray(array)
    band.SetNoDataValue(0)
    band.FlushCache()
    ds.FlushCache()
    ds = None


def test_canonicalize_hydgrpdcd_dual_policy() -> None:
    assert mod.canonicalize_hydgrpdcd("A", dual_group_policy="assume_d") == ("A", "ok")
    assert mod.canonicalize_hydgrpdcd("B/D", dual_group_policy="assume_d") == (
        "D",
        "ok_dual_assume_d",
    )
    assert mod.canonicalize_hydgrpdcd("C/D", dual_group_policy="error") == (
        None,
        "dual_group_blocked",
    )
    assert mod.canonicalize_hydgrpdcd(None, dual_group_policy="assume_d") == (
        None,
        "missing_raw",
    )
    assert mod.canonicalize_hydgrpdcd("W", dual_group_policy="assume_d") == (
        None,
        "unknown_raw_code",
    )


def test_build_lookup_records_status_and_codes() -> None:
    mukeys = [1, 2, 3, 4]
    raw = {1: "A", 2: "A/D", 3: None, 4: "W"}
    raw_source = {1: "muaggatt", 2: "component_hydgrp_major", 3: "muaggatt"}

    lut, records, status_counts, source_counts = mod.build_lookup_records(
        mukeys,
        raw,
        dual_group_policy="assume_d",
        raw_source_lookup=raw_source,
    )

    assert lut[1] == 1
    assert lut[2] == 4
    assert lut[3] == 0
    assert lut[4] == 0

    status_by_mukey = {row["mukey"]: row["status"] for row in records}
    assert status_by_mukey[1] == "ok"
    assert status_by_mukey[2] == "ok_dual_assume_d"
    assert status_by_mukey[3] == "missing_raw"
    assert status_by_mukey[4] == "unknown_raw_code"

    assert status_counts["ok"] == 1
    assert status_counts["ok_dual_assume_d"] == 1
    assert status_counts["missing_raw"] == 1
    assert status_counts["unknown_raw_code"] == 1
    assert source_counts["muaggatt"] == 2
    assert source_counts["component_hydgrp_major"] == 1
    assert source_counts["unknown"] == 1


def test_select_component_hydgrp_rows_prefers_first_non_empty() -> None:
    rows = [
        [100, 1, 70, None],
        [100, 2, 20, "B/D"],
        [100, 3, 10, "A"],
        [200, 4, 60, ""],
        [200, 5, 30, "C"],
        [300, 6, 80, None],
    ]
    selected = mod._select_component_hydgrp_rows(rows)
    assert selected[100] == "B/D"
    assert selected[200] == "C"
    assert selected[300] is None


def test_map_horizon_simple_texture_to_hsg() -> None:
    assert mod.map_horizon_simple_texture_to_hsg(sand_value=70, clay_value=10) == (
        "B",
        "texture4_mapped",
        "sand loam",
    )
    assert mod.map_horizon_simple_texture_to_hsg(sand_value=45, clay_value=15) == (
        "B",
        "texture4_mapped",
        "loam",
    )
    assert mod.map_horizon_simple_texture_to_hsg(sand_value=20, clay_value=20) == (
        "C",
        "texture4_mapped",
        "silt loam",
    )
    assert mod.map_horizon_simple_texture_to_hsg(sand_value=35, clay_value=35) == (
        "D",
        "texture4_mapped",
        "clay loam",
    )
    assert mod.map_horizon_simple_texture_to_hsg(sand_value=None, clay_value=20) == (
        None,
        "texture4_missing_sand_or_clay",
        "",
    )
    assert mod.map_horizon_simple_texture_to_hsg(sand_value=90, clay_value=20) == (
        None,
        "texture4_invalid_sand_clay",
        "",
    )


def test_select_chorizon_texture_fallback_rows_iterates() -> None:
    rows = [
        [10, 1001, 80, 5001, 0, 5, None, None, "VAR"],
        [10, 1001, 80, 5002, 5, 15, 70, 10, "SL"],
        [20, 2001, 70, 6001, 0, 10, 20, 20, "SIL"],
        [30, 3001, 90, 7001, 0, 10, None, None, None],
    ]
    selected = mod._select_chorizon_texture_fallback_rows(rows)
    assert selected[10][0] == "B"
    assert selected[20][0] == "C"
    assert selected[30][0] is None
    assert selected[30][1] == "texture4_missing_sand_or_clay"


def test_write_coded_raster_preserves_metadata_and_values(tmp_path: Path) -> None:
    source = tmp_path / "source.tif"
    output = tmp_path / "output.tif"
    source_data = np.array(
        [
            [0, 1, 2, 3],
            [4, 5, 2, 99],
        ],
        dtype=np.uint32,
    )
    _create_test_uint32_raster(source, source_data)

    lut = np.zeros(128, dtype=np.uint8)
    lut[1] = 1
    lut[2] = 2
    lut[3] = 3
    lut[4] = 4
    lut[5] = 4

    summary = mod.write_coded_raster(source, output, lut, nodata_code=0)
    assert summary["width"] == 4
    assert summary["height"] == 2
    assert summary["cell_counts_by_code"]["0"] == 2
    assert summary["cell_counts_by_code"]["1"] == 1
    assert summary["cell_counts_by_code"]["2"] == 2
    assert summary["cell_counts_by_code"]["3"] == 1
    assert summary["cell_counts_by_code"]["4"] == 2

    src_ds = gdal.Open(str(source), gdal.GA_ReadOnly)
    out_ds = gdal.Open(str(output), gdal.GA_ReadOnly)
    assert src_ds.GetGeoTransform() == out_ds.GetGeoTransform()
    assert src_ds.GetProjection() == out_ds.GetProjection()

    out_band = out_ds.GetRasterBand(1)
    assert out_band.DataType == gdal.GDT_Byte
    assert int(out_band.GetNoDataValue()) == 0

    expected = np.array(
        [
            [0, 1, 2, 3],
            [4, 4, 2, 0],
        ],
        dtype=np.uint8,
    )
    np.testing.assert_array_equal(out_band.ReadAsArray(), expected)


def test_fill_single_pixel_zero_holes_inplace(tmp_path: Path) -> None:
    raster = tmp_path / "holes.tif"
    source_data = np.array(
        [
            [2, 2, 2, 2, 2, 2],
            [2, 1, 1, 1, 2, 2],
            [2, 1, 0, 1, 2, 2],  # isolated hole -> should fill
            [2, 1, 1, 1, 0, 0],  # zero-cluster away from center -> should remain
            [2, 2, 2, 2, 0, 0],
        ],
        dtype=np.uint32,
    )
    _create_test_uint32_raster(raster, source_data)

    stats = mod.fill_single_pixel_zero_holes_inplace(raster)
    assert stats["filled_total"] == 1
    assert stats["filled_by_code"]["1"] == 1

    out_ds = gdal.Open(str(raster), gdal.GA_ReadOnly)
    out = out_ds.GetRasterBand(1).ReadAsArray()
    expected = np.array(
        [
            [2, 2, 2, 2, 2, 2],
            [2, 1, 1, 1, 2, 2],
            [2, 1, 1, 1, 2, 2],  # filled
            [2, 1, 1, 1, 0, 0],  # unchanged
            [2, 2, 2, 2, 0, 0],
        ],
        dtype=np.uint8,
    )
    np.testing.assert_array_equal(out, expected)


def test_apply_single_pixel_fill_to_summary() -> None:
    summary = {
        "width": 5,
        "height": 5,
        "cell_counts_by_code": {"0": 3, "1": 10, "2": 12, "3": 0, "4": 0},
    }
    fill = {
        "filled_total": 2,
        "filled_by_code": {"1": 1, "2": 1, "3": 0, "4": 0},
    }
    updated = mod.apply_single_pixel_fill_to_summary(summary, fill)
    assert updated["cell_counts_by_code"]["0"] == 1
    assert updated["cell_counts_by_code"]["1"] == 11
    assert updated["cell_counts_by_code"]["2"] == 13

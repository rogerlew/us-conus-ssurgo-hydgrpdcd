#!/usr/bin/env python3
"""Render a quicklook PNG map for the hydgrpdcd raster."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from osgeo import gdal


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render hydgrpdcd map PNG.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("/home/geodata/ssurgo/hydgrpdcd/hydgrpdcd.tif"),
        help="Input hydgrpdcd GeoTIFF path.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("/home/geodata/ssurgo/hydgrpdcd/hydgrpdcd_map.png"),
        help="Output PNG path.",
    )
    parser.add_argument(
        "--max-width",
        type=int,
        default=5000,
        help="Maximum raster width used for plotting.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="Output figure DPI.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    ds = gdal.Open(str(args.input))
    if ds is None:
        raise FileNotFoundError(f"Unable to open input raster: {args.input}")

    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    nodata = band.GetNoDataValue()
    if nodata is None:
        nodata = 0

    if args.max_width > 0 and arr.shape[1] > args.max_width:
        step = max(1, math.ceil(arr.shape[1] / args.max_width))
        arr = arr[::step, ::step]

    # Codebook: 0 nodata, 1 A, 2 B, 3 C, 4 D
    cmap = ListedColormap(
        [
            (1.0, 1.0, 1.0, 1.0),      # 0 NoData
            (0.60, 0.86, 0.55, 1.0),   # 1 A
            (0.98, 0.84, 0.50, 1.0),   # 2 B
            (0.96, 0.58, 0.35, 1.0),   # 3 C
            (0.80, 0.20, 0.20, 1.0),   # 4 D
        ]
    )

    counts = np.bincount(arr.ravel(), minlength=5).astype(np.int64)
    total_valid = max(1, int(np.sum(counts[1:5])))

    fig = plt.figure(figsize=(14, 9), dpi=args.dpi)
    ax = fig.add_subplot(1, 1, 1)
    ax.imshow(arr, cmap=cmap, vmin=0, vmax=4, interpolation="nearest")
    ax.set_axis_off()
    ax.set_title("CONUS SSURGO Hydrologic Soil Group (hydgrpdcd)", fontsize=15, weight="bold")

    legend_patches = [
        mpatches.Patch(color=cmap(1), label=f"A ({counts[1]:,}, {counts[1]/total_valid:.1%})"),
        mpatches.Patch(color=cmap(2), label=f"B ({counts[2]:,}, {counts[2]/total_valid:.1%})"),
        mpatches.Patch(color=cmap(3), label=f"C ({counts[3]:,}, {counts[3]/total_valid:.1%})"),
        mpatches.Patch(color=cmap(4), label=f"D ({counts[4]:,}, {counts[4]/total_valid:.1%})"),
    ]
    ax.legend(
        handles=legend_patches,
        loc="lower right",
        frameon=True,
        framealpha=0.95,
        fontsize=9,
        title="Hydrologic Group",
    )

    ax.text(
        0.01,
        0.01,
        f"NoData cells: {counts[0]:,}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[info] wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

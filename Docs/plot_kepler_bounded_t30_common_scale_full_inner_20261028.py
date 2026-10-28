#!/usr/bin/env python3
"""Show the t=30 Kepler density on one common scale, with an inner zoom."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle
import numpy as np

from plot_kepler_bounded_t0_t30_density_compare_20261028 import (
    CENTER,
    DOMAIN,
    FINAL_DEFAULT,
    INITIAL_DEFAULT,
    N_EXPECTED,
    raster,
    read_ordered,
    weighted_mean,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    here = Path(__file__).resolve().parent
    stem = "KEPLER_Q2_BOUNDED_W0_DENSITY_T30_COMMON_SCALE_FULL_INNER_2026-10-28"
    parser = argparse.ArgumentParser()
    parser.add_argument("--initial", type=Path, default=INITIAL_DEFAULT)
    parser.add_argument("--final", type=Path, default=FINAL_DEFAULT)
    parser.add_argument("--output", type=Path, default=here / f"{stem}.png")
    parser.add_argument("--summary", type=Path, default=here / f"{stem}.json")
    parser.add_argument("--pixels", type=int, default=1500)
    parser.add_argument("--percentile", type=float, default=97.5)
    args = parser.parse_args()
    if args.output.exists() or args.summary.exists():
        raise FileExistsError("Refusing to overwrite an existing output")
    if args.pixels < 100 or not 90.0 <= args.percentile < 100.0:
        raise ValueError("pixels must be >=100 and percentile in [90,100)")

    initial, meta0 = read_ordered(args.initial)
    final, meta1 = read_ordered(args.final)
    if len(initial) != N_EXPECTED or len(final) != N_EXPECTED:
        raise ValueError("Unexpected cell count")
    if abs(float(meta0["t"])) > 1e-12 or abs(float(meta1["t"]) - 30.0) > 1e-12:
        raise ValueError(f"Unexpected snapshot times: {meta0['t']} and {meta1['t']}")
    if initial["mass"].tobytes() != final["mass"].tobytes():
        raise ValueError("Cell masses are not bitwise equal by ID")
    if np.any(initial["w2"] != 0) or np.any(final["w2"] != 0):
        raise ValueError("Both endpoints must have W=0")

    xy0 = np.column_stack((initial["x"].astype(float), initial["y"].astype(float)))
    xy1 = np.column_stack((final["x"].astype(float), final["y"].astype(float)))
    rho0 = initial["den"].astype(float)
    rho1 = final["den"].astype(float)
    if (np.any(~np.isfinite(rho0)) or np.any(~np.isfinite(rho1))
            or np.any(rho0 <= 0.0) or np.any(rho1 <= 0.0)):
        raise ValueError("Density contains non-finite or non-positive values")

    r0 = np.hypot(xy0[:, 0] - CENTER[0], xy0[:, 1] - CENTER[1])
    disk = (r0 >= 0.5) & (r0 < 2.0)
    mass = initial["mass"].astype(float)
    rho_ref = weighted_mean(rho0[disk], mass[disk])
    contrast = 100.0 * (rho1 / rho_ref - 1.0)
    limit = float(np.quantile(np.abs(contrast[disk]), args.percentile / 100.0))
    if not np.isfinite(limit) or limit <= 0.0:
        raise ValueError("Density contrast has no positive range")

    full_extent = (*DOMAIN, *DOMAIN)
    zoom_half_width = 0.9
    zoom_extent = (CENTER[0] - zoom_half_width, CENTER[0] + zoom_half_width,
                   CENTER[1] - zoom_half_width, CENTER[1] + zoom_half_width)
    fields = []
    clipped = []
    for extent in (full_extent, zoom_extent):
        values, nearest = raster(xy1, contrast, extent, args.pixels)
        mask = disk[nearest]
        fields.append(np.ma.masked_where(~mask, values))
        clipped.append(int(np.count_nonzero(mask & (np.abs(values) > limit))))

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#171522")
    cmap.set_under("#164b9a")
    cmap.set_over("#a51224")
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 7.2), constrained_layout=True)
    titles = ("Full disk", "Inner-disk zoom")
    extents = (full_extent, zoom_extent)
    artists = []
    for ax, field, extent, title in zip(axes, fields, extents, titles):
        artists.append(ax.imshow(
            field, origin="lower", extent=extent, interpolation="nearest",
            cmap=cmap, norm=norm, rasterized=True,
        ))
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        radii = (0.5, 2.0) if extent == full_extent else (0.5, 0.75, 1.0)
        for radius in radii:
            ax.add_patch(Circle(CENTER, radius, fill=False, color="white",
                                linestyle="--", linewidth=0.75, alpha=0.85))

    fig.colorbar(
        artists[0], ax=axes.ravel().tolist(), orientation="horizontal",
        pad=0.10, shrink=0.76, aspect=48,
        label=rf"surface-density contrast $100(\Sigma/\Sigma_{{ref}}-1)$ [%]; "
              rf"identical scale in both panels, $\Sigma_{{ref}}={rho_ref:.8g}$",
        extend="both",
    )
    fig.suptitle(
        "Kepler surface density at t=30: shared scale and inner-disk detail\n"
        f"W=0, {N_EXPECTED:,} cells; ±{limit:.4g}% from the initial disk mean "
        f"({args.percentile:g}th percentile); nearest-generator raster, no smoothing",
        fontsize=14,
    )
    fig.savefig(args.output, dpi=200, facecolor="white")
    plt.close(fig)

    summary = {
        "method": "t=30 native W=0 nearest-generator density raster; contrast relative to initial mass-weighted disk mean",
        "initial": str(args.initial.resolve()),
        "final": str(args.final.resolve()),
        "initial_sha256": sha256(args.initial),
        "final_sha256": sha256(args.final),
        "initial_time": float(meta0["t"]),
        "final_time": float(meta1["t"]),
        "cell_count": N_EXPECTED,
        "fixed_mass_bitwise_equal_by_id": True,
        "W0_exact": True,
        "disk_initial_radius": [0.5, 2.0],
        "disk_cell_count": int(disk.sum()),
        "density_reference_initial_mass_weighted_disk_mean": rho_ref,
        "shared_symmetric_contrast_limit_percent": limit,
        "contrast_limit_percentile_by_cell": args.percentile,
        "cells_clipped_beyond_limit": int(np.count_nonzero(np.abs(contrast[disk]) > limit)),
        "cell_fraction_clipped": float(np.mean(np.abs(contrast[disk]) > limit)),
        "mass_weighted_mean_absolute_contrast_percent": weighted_mean(np.abs(contrast[disk]), mass[disk]),
        "raster_pixels_per_side": args.pixels,
        "pixel_samples_clipped_full_and_zoom": clipped,
        "inner_zoom_extent": list(zoom_extent),
        "rendering": "piecewise-constant nearest-generator assignment; identical linear norm and colorbar for full disk and zoom; no smoothing",
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"image={args.output.resolve()}")
    print(f"summary={args.summary.resolve()}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

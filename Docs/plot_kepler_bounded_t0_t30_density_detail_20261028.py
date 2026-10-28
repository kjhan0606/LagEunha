#!/usr/bin/env python3
"""Plot t=0 and t=30 Kepler density contrast on one scale, with an inner zoom."""

import argparse
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


def main():
    default_stem = "KEPLER_Q2_BOUNDED_W0_DENSITY_DETAIL_T0_VS_T30_2026-10-28"
    ap = argparse.ArgumentParser()
    ap.add_argument("--initial", type=Path, default=INITIAL_DEFAULT)
    ap.add_argument("--final", type=Path, default=FINAL_DEFAULT)
    ap.add_argument("--output", type=Path, default=Path(__file__).with_name(default_stem + ".png"))
    ap.add_argument("--summary", type=Path, default=Path(__file__).with_name(default_stem + ".json"))
    ap.add_argument("--pixels", type=int, default=1000)
    ap.add_argument("--percentile", type=float, default=99.5)
    args = ap.parse_args()
    if args.output.exists() or args.summary.exists():
        raise FileExistsError("Refusing to overwrite an existing output")
    if not 90.0 <= args.percentile < 100.0:
        raise ValueError("--percentile must be in [90, 100)")

    initial, meta0 = read_ordered(args.initial)
    final, meta1 = read_ordered(args.final)
    if len(initial) != N_EXPECTED or len(final) != N_EXPECTED:
        raise ValueError("Unexpected cell count")
    if abs(float(meta0["t"])) > 1e-12 or abs(float(meta1["t"]) - 30.0) > 1e-12:
        raise ValueError(f"Unexpected snapshot times: {meta0['t']} and {meta1['t']}")
    if initial["mass"].tobytes() != final["mass"].tobytes():
        raise ValueError("Cell masses are not bitwise equal by ID")
    if np.any(initial["w2"] != 0) or np.any(final["w2"] != 0):
        raise ValueError("Both snapshots must have W=0")

    xy0 = np.column_stack((initial["x"].astype(float), initial["y"].astype(float)))
    xy1 = np.column_stack((final["x"].astype(float), final["y"].astype(float)))
    rho0 = initial["den"].astype(float)
    rho1 = final["den"].astype(float)
    if np.any(~np.isfinite(rho0)) or np.any(~np.isfinite(rho1)) or np.any(rho0 <= 0) or np.any(rho1 <= 0):
        raise ValueError("Density contains non-finite or non-positive values")

    r0 = np.hypot(xy0[:, 0] - CENTER[0], xy0[:, 1] - CENTER[1])
    disk = (r0 >= 0.5) & (r0 < 2.0)
    mass = initial["mass"].astype(float)
    rho_ref = weighted_mean(rho0[disk], mass[disk])
    contrast0 = 100.0 * (rho0 / rho_ref - 1.0)
    contrast1 = 100.0 * (rho1 / rho_ref - 1.0)
    color_limit = float(np.quantile(
        np.abs(np.concatenate((contrast0[disk], contrast1[disk]))),
        args.percentile / 100.0,
    ))
    if color_limit <= 0:
        raise ValueError("Density contrast has no positive range")

    full_extent = (*DOMAIN, *DOMAIN)
    half_zoom = 1.0
    zoom_extent = (
        CENTER[0] - half_zoom, CENTER[0] + half_zoom,
        CENTER[1] - half_zoom, CENTER[1] + half_zoom,
    )

    def density_map(xy, contrast, extent):
        values, nearest = raster(xy, contrast, extent, args.pixels)
        mask = disk[nearest]
        return np.ma.masked_where(~mask, values)

    maps = (
        density_map(xy0, contrast0, full_extent),
        density_map(xy1, contrast1, full_extent),
        density_map(xy0, contrast0, zoom_extent),
        density_map(xy1, contrast1, zoom_extent),
    )
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("white")
    cmap.set_under("#183b8b")
    cmap.set_over("#8b1010")
    norm = TwoSlopeNorm(vmin=-color_limit, vcenter=0.0, vmax=color_limit)

    fig, axes = plt.subplots(2, 2, figsize=(12.8, 11.2), constrained_layout=True)
    extents = (full_extent, full_extent, zoom_extent, zoom_extent)
    titles = ("Initial density, full disk (t=0)", "Evolved density, full disk (t=30)",
              "Initial density, inner zoom (t=0)", "Evolved density, inner zoom (t=30)")
    images = []
    for ax, field, extent, title in zip(axes.flat, maps, extents, titles):
        images.append(ax.imshow(
            field, origin="lower", extent=extent, interpolation="nearest",
            cmap=cmap, norm=norm, rasterized=True,
        ))
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        for radius in ((0.5, 2.0) if extent == full_extent else (0.5, 1.0)):
            ax.add_patch(Circle(CENTER, radius, fill=False, color="black",
                                linestyle="--", linewidth=0.7, alpha=0.7))
    fig.colorbar(
        images[0], ax=axes.ravel().tolist(), pad=0.025, shrink=0.91, extend="both",
        label=rf"density contrast $100(\rho/\rho_{{ref}}-1)$ [%]; one shared scale, $\rho_{{ref}}={rho_ref:.8g}$",
    )
    fig.suptitle(
        f"Bounded aligned W=0 Kepler disk: density contrast at t=0 and t=30\n"
        f"{N_EXPECTED:,} cells; common symmetric color range ±{color_limit:.4g}% "
        f"({args.percentile:g}th percentile); nearest-generator raster, no smoothing",
        fontsize=14,
    )
    fig.savefig(args.output, dpi=180, facecolor="white")
    plt.close(fig)

    summary = {
        "method": "native W=0 nearest-generator raster; density contrast relative to initial mass-weighted disk mean",
        "initial": str(args.initial.resolve()),
        "final": str(args.final.resolve()),
        "initial_time": float(meta0["t"]),
        "final_time": float(meta1["t"]),
        "cell_count": N_EXPECTED,
        "fixed_mass_bitwise_equal_by_id": True,
        "W0_exact": True,
        "disk_initial_radius": [0.5, 2.0],
        "inner_zoom_extent": list(zoom_extent),
        "density_reference_mass_weighted_initial_disk_mean": rho_ref,
        "all_panels_shared_symmetric_density_contrast_limit_percent": color_limit,
        "limit_percentile_of_absolute_disk_contrast": args.percentile,
        "mass_weighted_mean_abs_density_contrast_percent_t0": weighted_mean(np.abs(contrast0[disk]), mass[disk]),
        "mass_weighted_mean_abs_density_contrast_percent_t30": weighted_mean(np.abs(contrast1[disk]), mass[disk]),
        "rendering": "piecewise-constant nearest-generator assignment; no smoothing; cells masked by initial disk cohort",
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"image={args.output.resolve()}")
    print(f"summary={args.summary.resolve()}")


if __name__ == "__main__":
    main()

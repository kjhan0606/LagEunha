#!/usr/bin/env python3
"""Make a same-scale t=0/t=30 density-field comparison focused on the inner disk."""

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
    FINAL_DEFAULT,
    INITIAL_DEFAULT,
    N_EXPECTED,
    read_ordered,
    raster,
    weighted_mean,
)


def main():
    stem = "KEPLER_Q2_BOUNDED_W0_DENSITY_T0_T30_INNER_FOCUS_2026-10-28"
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--initial", type=Path, default=INITIAL_DEFAULT)
    parser.add_argument("--final", type=Path, default=FINAL_DEFAULT)
    parser.add_argument("--output", type=Path, default=here / f"{stem}.png")
    parser.add_argument("--summary", type=Path, default=here / f"{stem}.json")
    parser.add_argument("--pixels", type=int, default=1400)
    parser.add_argument("--percentile", type=float, default=99.5)
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
        raise ValueError("Both snapshots must be W=0")

    xy0 = np.column_stack((initial["x"].astype(float), initial["y"].astype(float)))
    xy1 = np.column_stack((final["x"].astype(float), final["y"].astype(float)))
    rho0 = initial["den"].astype(float)
    rho1 = final["den"].astype(float)
    if (np.any(~np.isfinite(rho0)) or np.any(~np.isfinite(rho1))
            or np.any(rho0 <= 0.0) or np.any(rho1 <= 0.0)):
        raise ValueError("Density contains non-finite or non-positive values")

    r0 = np.hypot(xy0[:, 0] - CENTER[0], xy0[:, 1] - CENTER[1])
    disk = (r0 >= 0.5) & (r0 < 2.0)
    focus = (r0 >= 0.5) & (r0 < 1.0)
    inner_rim = (r0 >= 0.5) & (r0 < 0.75)
    mass = initial["mass"].astype(float)
    rho_ref = weighted_mean(rho0[disk], mass[disk])
    contrast0 = 100.0 * (rho0 / rho_ref - 1.0)
    contrast1 = 100.0 * (rho1 / rho_ref - 1.0)
    color_limit = float(np.quantile(
        np.abs(np.concatenate((contrast0[inner_rim], contrast1[inner_rim]))),
        args.percentile / 100.0,
    ))
    if not np.isfinite(color_limit) or color_limit <= 0.0:
        raise ValueError("Inner-disk density contrast has no positive range")

    half_width = 1.1
    extent = (CENTER[0] - half_width, CENTER[0] + half_width,
              CENTER[1] - half_width, CENTER[1] + half_width)

    def inner_map(xy, contrast):
        values, nearest = raster(xy, contrast, extent, args.pixels)
        return np.ma.masked_where(~focus[nearest], values)

    maps = (inner_map(xy0, contrast0), inner_map(xy1, contrast1))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("white")
    cmap.set_under("#183b8b")
    cmap.set_over("#8b1010")
    norm = TwoSlopeNorm(vmin=-color_limit, vcenter=0.0, vmax=color_limit)

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 6.8), constrained_layout=True)
    titles = ("Initial density field (t=0)", "Evolved density field (t=30)")
    images = []
    for ax, field, title in zip(axes, maps, titles):
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
        for radius in (0.5, 1.0):
            ax.add_patch(Circle(CENTER, radius, fill=False, color="#333333",
                                linestyle="--", linewidth=0.8, alpha=0.85))
    fig.colorbar(
        images[0], ax=axes.ravel().tolist(), pad=0.025, shrink=0.92, extend="both",
        label=rf"density contrast $100(\rho/\rho_{{ref}}-1)$ [%]; shared scale, "
              rf"$\rho_{{ref}}={rho_ref:.8g}$",
    )
    fig.suptitle(
        "Bounded aligned W=0 Kepler disk: inner density comparison\n"
        f"Same scale in both panels: ±{color_limit:.4g}% "
        f"(99.5th percentile of initial-radius 0.50–0.75 cells); "
        "shown cohort: initial radius 0.50–1.00; no smoothing",
        fontsize=13,
    )
    fig.savefig(args.output, dpi=200, facecolor="white")
    plt.close(fig)

    summary = {
        "method": "native W=0 nearest-generator raster of density contrast; no smoothing",
        "initial": str(args.initial.resolve()),
        "final": str(args.final.resolve()),
        "initial_time": float(meta0["t"]),
        "final_time": float(meta1["t"]),
        "cell_count": N_EXPECTED,
        "fixed_mass_bitwise_equal_by_id": True,
        "W0_exact": True,
        "reference_density": "initial mass-weighted mean for initial-radius 0.50<=r<2.00 cohort",
        "reference_density_value": rho_ref,
        "displayed_cohort_initial_radius": [0.5, 1.0],
        "color_limit_cohort_initial_radius": [0.5, 0.75],
        "shared_symmetric_color_limit_percent": color_limit,
        "color_limit_percentile": args.percentile,
        "inner_rim_cell_count": int(inner_rim.sum()),
        "inner_rim_mass_weighted_abs_contrast_percent_t0": weighted_mean(
            np.abs(contrast0[inner_rim]), mass[inner_rim]),
        "inner_rim_mass_weighted_abs_contrast_percent_t30": weighted_mean(
            np.abs(contrast1[inner_rim]), mass[inner_rim]),
        "rendering": "piecewise-constant nearest-generator assignment; pixel colors use stored cell density",
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"image={args.output.resolve()}")
    print(f"summary={args.summary.resolve()}")


if __name__ == "__main__":
    main()

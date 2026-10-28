#!/usr/bin/env python3
"""Compare t=0 and t=30 density fields for the bounded aligned W=0 Kepler run.

Each density raster assigns the nearest generator's native density to a pixel,
which visualizes the W=0 Voronoi cells without spatial smoothing. The t=0 and
t=30 panels share one linear color normalization. A third panel shows the
ID-matched fractional density change at the t=30 cell positions.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Circle
import numpy as np
from scipy.spatial import cKDTree


SNAPSHOT_MODULE = Path(
    "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020"
)
INITIAL_DEFAULT = Path(
    "/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/"
    "t001_v2/khout.000000.dat"
)
FINAL_DEFAULT = Path(
    "/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/"
    "wallselect_to30_a10_v1/khout.005247.dat"
)
N_EXPECTED = 55351
MASK56 = np.uint64((1 << 56) - 1)
DOMAIN = (0.75, 5.25)
CENTER = (3.0, 3.0)

sys.path.insert(0, str(SNAPSHOT_MODULE))
from audit_snapshot import read_snapshot  # noqa: E402


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_ordered(path: Path):
    data, meta = read_snapshot(path)
    if len(data) != N_EXPECTED:
        raise ValueError(f"Expected {N_EXPECTED} cells in {path}, found {len(data)}")
    ids = (data["u4if"].astype(np.uint64) & MASK56).astype(np.int64)
    order = np.argsort(ids)
    if not np.array_equal(ids[order], np.arange(N_EXPECTED)):
        raise ValueError(f"Cell IDs are incomplete or duplicated in {path}")
    return np.asarray(data[order]), meta


def raster(xy: np.ndarray, values: np.ndarray, extent, pixels: int):
    x0, x1, y0, y1 = extent
    x = np.linspace(x0, x1, pixels)
    y = np.linspace(y0, y1, pixels)
    gx, gy = np.meshgrid(x, y)
    nearest = cKDTree(xy).query(
        np.column_stack((gx.ravel(), gy.ravel())), k=1
    )[1]
    return np.asarray(values)[nearest].reshape(pixels, pixels), nearest.reshape(pixels, pixels)


def add_reference_radii(ax, radii, color="white"):
    for radius in radii:
        ax.add_patch(
            Circle(CENTER, radius, fill=False, color=color,
                   linestyle="--", linewidth=0.75, alpha=0.8)
        )


def weighted_mean(values, weights):
    return float(
        np.sum(values * weights, dtype=np.longdouble)
        / np.sum(weights, dtype=np.longdouble)
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--initial", type=Path, default=INITIAL_DEFAULT)
    ap.add_argument("--final", type=Path, default=FINAL_DEFAULT)
    ap.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name(
            "KEPLER_Q2_BOUNDED_W0_DENSITY_T0_VS_T30_2026-10-28.png"
        ),
    )
    ap.add_argument(
        "--summary",
        type=Path,
        default=Path(__file__).with_name(
            "KEPLER_Q2_BOUNDED_W0_DENSITY_T0_VS_T30_2026-10-28.json"
        ),
    )
    ap.add_argument("--pixels", type=int, default=1100)
    args = ap.parse_args()
    if args.output.exists() or args.summary.exists():
        raise FileExistsError("Refusing to overwrite an existing output")

    initial, meta0 = read_ordered(args.initial)
    final, meta1 = read_ordered(args.final)
    if abs(float(meta0["t"])) > 1e-12 or abs(float(meta1["t"]) - 30.0) > 1e-12:
        raise ValueError(f"Unexpected snapshot times: {meta0['t']} and {meta1['t']}")
    if initial["mass"].tobytes() != final["mass"].tobytes():
        raise ValueError("Cell masses are not bitwise equal by ID")
    if np.any(initial["w2"] != 0) or np.any(final["w2"] != 0):
        raise ValueError("Both endpoints must be W=0")

    xy0 = np.column_stack((initial["x"].astype(float), initial["y"].astype(float)))
    xy1 = np.column_stack((final["x"].astype(float), final["y"].astype(float)))
    rho0 = initial["den"].astype(float)
    rho1 = final["den"].astype(float)
    if not (np.all(np.isfinite(rho0)) and np.all(np.isfinite(rho1))):
        raise ValueError("Density contains non-finite values")
    if np.any(rho0 <= 0) or np.any(rho1 <= 0):
        raise ValueError("Density contains non-positive values")

    r0 = np.hypot(xy0[:, 0] - CENTER[0], xy0[:, 1] - CENTER[1])
    disk = (r0 >= 0.5) & (r0 < 2.0)
    delta = rho1 / rho0 - 1.0
    mass = initial["mass"].astype(float)
    density_limits = np.quantile(
        np.concatenate((rho0[disk], rho1[disk])), [0.01, 0.99]
    )
    change_limit = float(np.quantile(np.abs(delta[disk]), 0.99))
    if change_limit <= 0:
        raise ValueError("The ID-matched disk density change is identically zero")

    summary = {
        "method": "native W=0 nearest-generator raster; ID-matched Lagrangian density change",
        "initial": str(args.initial.resolve()),
        "final": str(args.final.resolve()),
        "initial_sha256": digest(args.initial),
        "final_sha256": digest(args.final),
        "initial_time": float(meta0["t"]),
        "final_time": float(meta1["t"]),
        "cell_count": N_EXPECTED,
        "fixed_mass_bitwise_equal_by_id": True,
        "W0_exact": True,
        "disk_initial_radius": [0.5, 2.0],
        "disk_cell_count": int(disk.sum()),
        "shared_density_color_limits_1st_99th_percentile": density_limits.tolist(),
        "density_change_symmetric_color_limit_99th_percentile": change_limit,
        "mass_weighted_abs_fractional_density_change": weighted_mean(np.abs(delta[disk]), mass[disk]),
        "mass_weighted_signed_fractional_density_change": weighted_mean(delta[disk], mass[disk]),
        "radial_bins": [],
        "rendering": "piecewise-constant nearest-generator assignment; no smoothing",
    }
    for label, lo, hi in (
        ("inner_0p50_0p65", 0.50, 0.65),
        ("inner_0p65_0p75", 0.65, 0.75),
        ("disk_mid_0p75_1p75", 0.75, 1.75),
        ("disk_outer_1p75_2p00", 1.75, 2.00),
    ):
        sel = (r0 >= lo) & (r0 < hi)
        summary["radial_bins"].append({
            "name": label,
            "initial_radius": [lo, hi],
            "cell_count": int(sel.sum()),
            "mass_weighted_abs_fractional_density_change": weighted_mean(np.abs(delta[sel]), mass[sel]),
            "mass_weighted_signed_fractional_density_change": weighted_mean(delta[sel], mass[sel]),
            "p50_abs_fractional_density_change": float(np.quantile(np.abs(delta[sel]), 0.50)),
            "p99_abs_fractional_density_change": float(np.quantile(np.abs(delta[sel]), 0.99)),
            "max_abs_fractional_density_change": float(np.max(np.abs(delta[sel]))),
        })

    full_extent = (*DOMAIN, *DOMAIN)
    map0, _ = raster(xy0, rho0, full_extent, args.pixels)
    map1, nearest1 = raster(xy1, rho1, full_extent, args.pixels)
    delta_map = delta[nearest1]
    disk_map = disk[nearest1]
    delta_map = np.ma.masked_where(~disk_map, delta_map)
    density_norm = Normalize(*density_limits, clip=False)
    density_cmap = plt.get_cmap("viridis").copy()
    density_cmap.set_under("#171522")
    density_cmap.set_over("#fff36b")
    change_norm = TwoSlopeNorm(vmin=-change_limit, vcenter=0.0, vmax=change_limit)

    fig, axes = plt.subplots(1, 3, figsize=(16.4, 5.9), constrained_layout=True)
    im0 = axes[0].imshow(
        map0, origin="lower", extent=full_extent, interpolation="nearest",
        cmap=density_cmap, norm=density_norm, rasterized=True,
    )
    axes[0].set_title(r"Initial density $t=0$")
    im1 = axes[1].imshow(
        map1, origin="lower", extent=full_extent, interpolation="nearest",
        cmap=density_cmap, norm=density_norm, rasterized=True,
    )
    axes[1].set_title(r"Evolved density $t=30$")
    im2 = axes[2].imshow(
        delta_map, origin="lower", extent=full_extent, interpolation="nearest",
        cmap="coolwarm", norm=change_norm, rasterized=True,
    )
    axes[2].set_title(r"ID-matched $\Sigma_{30}/\Sigma_0-1$")
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlim(DOMAIN)
        ax.set_ylim(DOMAIN)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        add_reference_radii(ax, (0.5, 2.0))
    fig.colorbar(
        im1, ax=axes[:2], pad=0.02, shrink=0.88, extend="both",
        label=r"surface density $\Sigma$ (shared linear scale for $t=0$ and $t=30$)",
    )
    fig.colorbar(
        im2, ax=axes[2], pad=0.02, shrink=0.88, extend="both",
        label=r"fractional density change (disk cells; symmetric 99th-percentile scale)",
    )
    fig.suptitle(
        r"Bounded aligned W=0 Kepler disk: density texture versus evolution" +
        f"\n{N_EXPECTED:,} cells; density endpoints use one shared 1st–99th percentile scale",
        fontsize=14,
    )
    fig.savefig(args.output, dpi=180, facecolor="white")
    plt.close(fig)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")

    print(f"image={args.output.resolve()}")
    print(f"summary={args.summary.resolve()}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

"""Render native polar Kepler density and exact-orbit Voronoi control."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.spatial import cKDTree

from audit_kepler_polar_checkpoint_20261026 import (
    CX, CY, disk_ballistic_areas, load_ordered,
)

DOMAIN = (.75, 5.25)
PIXELS = 720


def sample(xy, field, grid):
    nearest = cKDTree(xy).query(grid, k=1)[1]
    return np.asarray(field)[nearest].reshape(PIXELS, PIXELS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("initial", type=Path)
    ap.add_argument("endpoint", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    p0, m0 = load_ordered(args.initial)
    p, meta = load_ordered(args.endpoint)
    t = float(meta["t"])
    if m0["t"] != 0 or not (0 < t <= .825 + 1e-12):
        raise ValueError("Expected native t=0 and early endpoint")
    xy0 = np.column_stack((p0["x"], p0["y"]))
    xy_native = np.column_stack((p["x"], p["y"]))
    r0 = np.hypot(p0["x"]-CX, p0["y"]-CY)
    disk = (r0 >= .5) & (r0 < 2.)
    selected = np.flatnonzero(disk)
    xy_exact, area_exact = disk_ballistic_areas(p0, t, selected)
    native_change = np.full(len(p), np.nan)
    native_change[disk] = abs(p["den"][disk].astype(float) /
                              p0["den"][disk].astype(float) - 1)
    geom_change = np.full(len(p), np.nan)
    geom_change[selected] = abs(p0["volume"][selected].astype(float) /
                                area_exact - 1)
    axis = np.linspace(*DOMAIN, PIXELS)
    gx, gy = np.meshgrid(axis, axis)
    grid = np.column_stack((gx.ravel(), gy.ravel()))
    initial_density = sample(xy0, p0["den"], grid)
    native_density = sample(xy_native, p["den"], grid)
    native_error = sample(xy_native, native_change, grid)
    geometry_error = sample(xy_exact, geom_change, grid)

    fig, ax = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)
    for col, (title, density) in enumerate((("Native initial density", initial_density),
                                            (f"Native density at t={t:.6g}", native_density))):
        artist = ax[0, col].imshow(density, origin="lower", extent=(*DOMAIN, *DOMAIN),
                                   norm=LogNorm(vmin=.01, vmax=1.1, clip=True),
                                   cmap="viridis", interpolation="nearest")
        ax[0, col].set_title(title)
    fig.colorbar(artist, ax=ax[0, :], label="surface density")
    for col, (title, error) in enumerate((("Native |density / initial - 1|", native_error),
                                          ("Exact-orbit Voronoi prediction", geometry_error))):
        artist = ax[1, col].imshow(np.ma.masked_invalid(np.maximum(error, 1e-7)),
                                   origin="lower", extent=(*DOMAIN, *DOMAIN),
                                   norm=LogNorm(vmin=1e-6, vmax=1e-1, clip=True),
                                   cmap="magma", interpolation="nearest")
        ax[1, col].set_title(title)
    fig.colorbar(artist, ax=ax[1, :], label="absolute relative density change; disk cells")
    for panel in ax.flat:
        panel.set_aspect("equal")
        panel.set_xlabel("x")
        panel.set_ylabel("y")
        panel.set_xticks([1, 2, 3, 4, 5])
        panel.set_yticks([1, 2, 3, 4, 5])
    fig.suptitle("384² Kepler polar IC: native RK4+CUDA versus orbit-only geometry\n"
                 "Early checkpoint; not the GIZMO long-time endpoint")
    fig.savefig(args.output, dpi=180)
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()

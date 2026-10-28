"""Render native W=0 Kepler density at any complete checkpoint."""

import argparse
from pathlib import Path
import struct
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot

CX = CY = 3.
DOMAIN = (.75, 5.25)
PIXELS = 720


def complete_ordered(path):
    with path.open("rb") as stream:
        magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
    if (magic, stride, stress) != (-1, 576, 1):
        raise ValueError(f"Unexpected native snapshot header: {path}")
    if path.stat().st_size != 16 + n*stride + 52:
        raise ValueError(f"Incomplete native snapshot: {path}")
    p, meta = read_snapshot(path)
    ids = (p["u4if"].astype(np.uint64) & np.uint64((1 << 56)-1)).astype(np.int64)
    if not np.array_equal(np.sort(ids), np.arange(n)):
        raise ValueError(f"Incomplete cell IDs: {path}")
    return p[np.argsort(ids)], meta


def raster(xy, field, grid):
    ids = cKDTree(xy).query(grid, k=1, workers=1)[1]
    return np.asarray(field)[ids].reshape(PIXELS, PIXELS)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("endpoint", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    p0, m0 = complete_ordered(args.initial)
    p, mt = complete_ordered(args.endpoint)
    if len(p0) != len(p) or float(m0["t"]) != 0 or float(mt["t"]) <= 0:
        raise ValueError("Expected matching native t=0 and later snapshots")
    if not np.array_equal(p0["mass"], p["mass"]):
        raise ValueError("Cell masses changed")
    if not (np.all(p0["w2"] == 0) and np.all(p["w2"] == 0)):
        raise ValueError("Nearest-site raster requires W=0 Voronoi geometry")
    if not np.all(np.isfinite(p0["den"]) & (p0["den"] > 0) &
                  np.isfinite(p["den"]) & (p["den"] > 0)):
        raise ValueError("Invalid density")
    xy0 = np.column_stack((p0["x"], p0["y"]))
    xy = np.column_stack((p["x"], p["y"]))
    r0 = np.hypot(p0["x"]-CX, p0["y"]-CY)
    disk = (r0 >= .5) & (r0 < 2.)
    density_ratio = p["den"].astype(float)/p0["den"].astype(float)-1.
    density_ratio[~disk] = np.nan
    axis = np.linspace(*DOMAIN, PIXELS)
    gx, gy = np.meshgrid(axis, axis)
    grid = np.column_stack((gx.ravel(), gy.ravel()))
    image0 = raster(xy0, p0["den"], grid)
    image = raster(xy, p["den"], grid)
    change = raster(xy, density_ratio, grid)

    fig, panels = plt.subplots(1, 3, figsize=(15, 4.7), constrained_layout=True)
    for ax, field, title in ((panels[0], image0, "Initial native density"),
                             (panels[1], image, f"Native density, t={float(mt['t']):.6g}")):
        artist = ax.imshow(field, origin="lower", extent=(*DOMAIN, *DOMAIN),
                           norm=LogNorm(vmin=.01, vmax=1.1, clip=True),
                           cmap="viridis", interpolation="nearest")
        ax.set_title(title)
    fig.colorbar(artist, ax=panels[:2], label="surface density", shrink=.85)
    artist = panels[2].imshow(np.ma.masked_invalid(change), origin="lower",
                              extent=(*DOMAIN, *DOMAIN),
                              norm=SymLogNorm(linthresh=1e-4, linscale=1.,
                                              vmin=-.1, vmax=.1, base=10),
                              cmap="coolwarm", interpolation="nearest")
    panels[2].set_title("Disk density / initial - 1")
    fig.colorbar(artist, ax=panels[2], label="relative change", shrink=.85)
    for ax in panels:
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.suptitle(f"Native CUDA+MPI+OpenMP Kepler, {len(p):,} fixed-mass W=0 cells")
    fig.savefig(args.output, dpi=180)
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()

"""Render the accepted C010 t=0 and t=6 W=0 Voronoi density fields."""

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/lagforce_matplotlib_20261026")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path("/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
sys.path.insert(0, str(ROOT))
from audit_snapshot import read_snapshot

SOURCES = [
    ROOT / "core_refinement/fine_integrated/khout.000001.dat",
    ROOT / "rk4_diagnostics/varpress_w0_geom_cfl_long_20261026/C010/to6/khout.032479.dat",
]
OUTPUT = Path(__file__).with_name("KEPLER_C010_DENSITY_T0_T6_2026-10-26.png")

grid = np.linspace(3, 21, 600)
xx, yy = np.meshgrid(grid, grid, indexing="xy")
query = np.column_stack((xx.ravel(), yy.ravel()))

fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)
for ax, source in zip(axes, SOURCES):
    p, meta = read_snapshot(source)
    seeds = np.column_stack((p["x"], p["y"]))
    nearest = cKDTree(seeds).query(query, workers=1)[1]
    rho = p["den"][nearest].reshape(xx.shape)
    im = ax.imshow(rho, origin="lower", extent=(3, 21, 3, 21),
                   cmap="viridis", norm=LogNorm(vmin=1e-3, vmax=1.2),
                   interpolation="nearest", rasterized=True)
    t = float(meta["t"])
    ax.set_title(f"C010 density at t={t:g} ({t/17.778891990936625:.3f} inner orbits)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")

fig.colorbar(im, ax=axes, label="Cell density", shrink=0.78)
fig.suptitle("Historical tapered 1/r disk, 71,080 W=0 cells; not the GIZMO initial profile")
fig.savefig(OUTPUT, dpi=180)
print(OUTPUT)

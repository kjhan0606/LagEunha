"""Render real power-cell density fields for the accepted t=8.2507 snapshots."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.path import Path as MplPath
import numpy as np

from plot_kepler_q2_weighted_vs_w0_density_20261027 import DOMAIN, case_fields
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_weighted_t10_a10_v1")
RECEIPT = DOCS / "KEPLER_Q2_WEIGHTED_STEP000903_MONITOR_2026-10-27.json"
OUTPUT = DOCS / "KEPLER_Q2_WEIGHTED_ALIGNED_JITTER_DENSITY_T825_2026-10-27.png"


def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    receipt = json.loads(RECEIPT.read_text())
    if receipt["step"] != 903:
        raise ValueError("Wrong checkpoint receipt")
    axis = np.linspace(*DOMAIN, 512)
    xx, yy = np.meshgrid(axis, axis)
    points = np.column_stack((xx.ravel(), yy.ravel()))
    theta = np.arange(1024)*2*np.pi/1024
    polygon = np.column_stack((3+2.75*np.cos(theta), 3+2.75*np.sin(theta)))
    inside = MplPath(polygon).contains_points(points).reshape(xx.shape)
    fig, panels = plt.subplots(2, 2, figsize=(10.5, 9), constrained_layout=True)
    times = []
    for column, label in enumerate(("aligned", "jitter")):
        run = ROOT / label
        start, end = run / "initial.dat", run / "khout.000903.dat"
        row = receipt["cases"][label]
        if (digest(start) != row["initial_sha256"]
                or digest(end) != row["weighted_endpoint_sha256"]):
            raise ValueError(f"{label}: checkpoint changed")
        t = snapshot_time(end, 55351)
        if abs(t-row["time"]) > 1e-10:
            raise ValueError(f"{label}: wrong checkpoint time")
        times.append(t)
        rho, relative, disk_abs = case_fields(start, end, points, weighted=True)
        if abs(disk_abs-row["bands"]["disk"]["weighted"]["density_abs_mass_mean"]) > 1e-10:
            raise ValueError(f"{label}: density metric mismatch")
        density = np.ma.masked_where(~inside, rho.reshape(xx.shape))
        residual = np.ma.masked_invalid(np.where(
            inside, relative.reshape(xx.shape), np.nan))
        artist_density = panels[0, column].imshow(
            density, origin="lower", extent=(*DOMAIN, *DOMAIN),
            norm=LogNorm(vmin=.01, vmax=1.1, clip=True), cmap="viridis",
            interpolation="nearest")
        artist_residual = panels[1, column].imshow(
            residual, origin="lower", extent=(*DOMAIN, *DOMAIN),
            norm=SymLogNorm(linthresh=1e-5, vmin=-.01, vmax=.01, base=10),
            cmap="coolwarm", interpolation="nearest")
        panels[0, column].set_title(f"{label}: surface density")
        panels[1, column].set_title(f"disk mean |ΔΣ/Σ₀| = {disk_abs:.4g}")
    if abs(times[0]-times[1]) > 1e-9:
        raise ValueError("Unequal IC snapshot times")
    fig.colorbar(artist_density, ax=panels[0, :].tolist(), label="surface density Σ")
    fig.colorbar(artist_residual, ax=panels[1, :].tolist(), label="Σ/Σ₀ − 1")
    for panel in panels.flat:
        panel.set_aspect("equal")
        panel.set_xlabel("x")
        panel.set_ylabel("y")
        panel.set_xticks((1, 2, 3, 4, 5))
        panel.set_yticks((1, 2, 3, 4, 5))
    fig.suptitle(f"Native bounded q2 weighted Kepler density, t={times[0]:.6f}")
    fig.savefig(OUTPUT, dpi=150)
    plt.close(fig)
    print(OUTPUT)


if __name__ == "__main__":
    main()

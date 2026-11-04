"""Plot actual periodic W=0 Voronoi density regions in completed native KH runs."""
from pathlib import Path
import argparse
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

from rotating_shear_assessment import read_native

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("suite", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
args.output.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.labelsize": 10,
                     "xtick.labelsize": 9, "ytick.labelsize": 9})
fig, axes = plt.subplots(1, 3, figsize=(7.2, 4.8), sharex=True, sharey=True,
                         gridspec_kw={"wspace": 0})
x = (np.arange(512) + 0.5) / 512
y = 2 * (np.arange(1024) + 0.5) / 1024
xx, yy = np.meshgrid(x, y)
samples = np.column_stack([xx.ravel(), yy.ravel()])
metadata = []
for ax, name, label in zip(axes,
        ["kh_64x128", "kh_128x256", "kh_256x512"],
        ["64 × 128", "128 × 256", "256 × 512"]):
    directory = args.suite / name
    result = json.loads((directory / "result.json").read_text())
    assert result["exit"] == 0 and result["latest_time"] == 2
    path = next(Path(row["path"]) for row in result["outputs"] if row["time"] == 2)
    state, time = read_native(path)
    assert np.all(state["W"] == 0) and np.all(state["pi"] == 0)
    points = np.column_stack([np.mod(state["x"], 1), np.mod(state["y"], 2)])
    rho = state["m"] / state["V"]
    tree = cKDTree(points, boxsize=[1, 2])
    _, nearest = tree.query(samples, workers=4)
    image = ax.imshow(rho[nearest].reshape(xx.shape), origin="lower",
                      extent=[0, 1, 0, 2], cmap="viridis", vmin=0.9, vmax=2.2,
                      interpolation="nearest", aspect="equal")
    ax.text(0.04, 1.94, label, va="top", ha="left", color="white", fontsize=10,
            bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none", "pad": 2})
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 2)
    ax.set_xlabel("x")
    ax.set_xticks([0, 0.5, 1])
    ax.tick_params(direction="in", top=True, right=True)
    metadata.append({"case": name, "native_snapshot": str(path), "time": time,
                     "density_min": float(rho.min()), "density_max": float(rho.max()),
                     "cells_outside_colour_limits": int(np.sum((rho < 0.9) | (rho > 2.2)))})
# Retain every panel's axes and tick marks, suppress only overlapping labels.
for i, ax in enumerate(axes):
    labels = ax.get_xticklabels()
    if i:
        labels[0].set_visible(False)
    if i < len(axes) - 1:
        labels[-1].set_visible(False)
axes[0].set_ylabel("y")
fig.subplots_adjust(left=0.08, right=0.87, bottom=0.12, top=0.98, wspace=0)
bar = fig.colorbar(image, cax=fig.add_axes([0.9, 0.16, 0.025, 0.78]), ticks=[1, 1.5, 2])
bar.set_label(r"$\rho$")
fig.savefig(args.output / "KH_density_t2.png", dpi=200)
fig.savefig(args.output / "KH_density_t2.pdf")
plt.close(fig)
(args.output / "metadata.json").write_text(json.dumps({
    "scope": "Native common-operator periodic KH at t=2. Pixels sample exact nearest-generator W0 Voronoi regions, rho=m/V. No interpolation between cells; no independent-reference accuracy claim.",
    "colour_limits": [0.9, 2.2], "cases": metadata}, indent=2) + "\n")

"""Compare accepted bounded Kepler density fields with common color scales."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.path import Path as MplPath
import numpy as np
from scipy.spatial import cKDTree

from plot_kepler_bounded_periodic_density_20261027 import (
    DOMAIN, PIXELS, POLYGON_AREA, load_case,
)
from prepare_kepler_sparse_area_restart_20261027 import digest


def accepted_case(label, receipt_path, expected_failed=()):
    receipt = json.loads(receipt_path.read_text())
    failed = {key for key, passed in receipt["checks"].items() if not passed}
    if failed != set(expected_failed):
        raise ValueError(f"{label}: native failed checks {sorted(failed)} differ from expected {sorted(set(expected_failed))}")
    initial = Path(receipt["initial"])
    endpoint = Path(receipt["endpoint"])
    if digest(initial) != receipt["initial_sha256"] or digest(endpoint) != receipt["endpoint_sha256"]:
        raise ValueError(f"{label}: checkpoint differs from accepted receipt")
    case = load_case(label, initial, endpoint, POLYGON_AREA)
    if abs(case["time"] - float(receipt["time"])) > 1e-10:
        raise ValueError(f"{label}: time differs from accepted receipt")
    if abs(case["disk_abs"] - receipt["disk_density_abs_mass_weighted"]) > 1e-10:
        raise ValueError(f"{label}: disk density differs from accepted receipt")
    return case


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("aligned_receipt", type=Path)
    parser.add_argument("jitter_receipt", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jitter-expected-failed-check", action="append", default=[])
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    cases = (
        accepted_case("Aligned rings", args.aligned_receipt),
        accepted_case("5% independent site jitter", args.jitter_receipt,
                      args.jitter_expected_failed_check),
    )
    if abs(cases[0]["time"] - cases[1]["time"]) > 1e-10:
        raise ValueError("The two accepted states have different physical times")
    axis = np.linspace(*DOMAIN, PIXELS)
    xx, yy = np.meshgrid(axis, axis)
    query = np.column_stack((xx.ravel(), yy.ravel()))
    theta = np.arange(1024) * 2 * np.pi / 1024
    boundary = np.column_stack((3 + 2.75 * np.cos(theta),
                                3 + 2.75 * np.sin(theta)))
    inside = MplPath(boundary).contains_points(query).reshape(PIXELS, PIXELS)
    fig, panels = plt.subplots(2, 2, figsize=(10.8, 9.2), constrained_layout=True)
    for col, case in enumerate(cases):
        nearest = cKDTree(case["xy"]).query(query, k=1, workers=1)[1]
        density = case["density"][nearest].reshape(PIXELS, PIXELS)
        relative = case["relative"][nearest].reshape(PIXELS, PIXELS)
        density = np.ma.masked_where(~inside, density)
        relative = np.ma.masked_where(~inside, relative)
        top = panels[0, col].imshow(
            density, origin="lower", extent=(*DOMAIN, *DOMAIN),
            norm=LogNorm(vmin=.01, vmax=1.1, clip=True),
            cmap="viridis", interpolation="nearest")
        bottom = panels[1, col].imshow(
            relative, origin="lower", extent=(*DOMAIN, *DOMAIN),
            norm=SymLogNorm(linthresh=1e-5, linscale=1., vmin=-.005,
                            vmax=.005, base=10), cmap="coolwarm",
            interpolation="nearest")
        panels[0, col].set_title(case["label"])
        panels[1, col].set_title(f"disk mean |ΔΣ/Σ₀| = {case['disk_abs']:.3g}")
    fig.colorbar(top, ax=panels[0, :].tolist(), label="surface density Σ")
    fig.colorbar(bottom, ax=panels[1, :].tolist(), label="disk Σ/Σ₀ − 1")
    for panel in panels.flat:
        panel.set_aspect("equal")
        panel.set_xlabel("x")
        panel.set_ylabel("y")
        panel.set_xticks((1, 2, 3, 4, 5))
        panel.set_yticks((1, 2, 3, 4, 5))
    failure_label = ("; exploratory, jitter FAILED "
                     + ", ".join(sorted(args.jitter_expected_failed_check))
                     if args.jitter_expected_failed_check else "")
    fig.suptitle(f"Native bounded Kepler density at t={cases[0]['time']:.6g}; "
                 f"common [-2,2] view and colors{failure_label}")
    fig.savefig(args.output, dpi=150)
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()

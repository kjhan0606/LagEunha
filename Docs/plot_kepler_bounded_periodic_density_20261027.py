"""Show native q2 periodic and fixed-polygon density on common paper scales."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.path import Path as MplPath
import numpy as np
from scipy.spatial import cKDTree

from plot_kepler_native_density_long_20261026 import complete_ordered


DOMAIN = (1.0, 5.0)
PIXELS = 720
POLYGON_AREA = 0.5 * 1024 * 2.75**2 * np.sin(2 * np.pi / 1024)


def load_case(label, initial, endpoint, expected_area):
    p0, m0 = complete_ordered(initial)
    p1, m1 = complete_ordered(endpoint)
    if len(p0) != len(p1) or abs(float(m0["t"])) > 1e-12:
        raise ValueError(f"{label}: invalid initial/native endpoint")
    if not np.array_equal(p0["mass"], p1["mass"]):
        raise ValueError(f"{label}: mass changed")
    if np.any(p0["w2"] != 0) or np.any(p1["w2"] != 0):
        raise ValueError(f"{label}: requires W=0 nearest-site cells")
    for p in (p0, p1):
        if np.any(~np.isfinite(p["den"])) or np.any(p["den"] <= 0):
            raise ValueError(f"{label}: invalid density")
        if abs(np.sum(p["volume"].astype(float), dtype=np.longdouble)
               - expected_area) > 1e-4:
            raise ValueError(f"{label}: cell areas do not tile their domain")
    r0 = np.hypot(p0["x"].astype(float) - 3, p0["y"].astype(float) - 3)
    disk = (r0 >= 0.5) & (r0 < 2.0)
    relative = p1["den"].astype(float) / p0["den"].astype(float) - 1
    mass = p0["mass"][disk].astype(float)
    disk_abs = float(np.sum(abs(relative[disk]) * mass, dtype=np.longdouble)
                     / np.sum(mass, dtype=np.longdouble))
    relative[~disk] = np.nan
    return {"label": label, "time": float(m1["t"]), "count": len(p1),
            "xy": np.column_stack((p1["x"], p1["y"])),
            "density": p1["den"].astype(float), "relative": relative,
            "disk_abs": disk_abs}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("periodic_initial", type=Path)
    parser.add_argument("periodic_endpoint", type=Path)
    parser.add_argument("bounded_initial", type=Path)
    parser.add_argument("bounded_endpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    cases = (
        load_case("Sparse q2 periodic", args.periodic_initial,
                  args.periodic_endpoint, 36.0),
        load_case("Sparse q2 fixed 1024-gon", args.bounded_initial,
                  args.bounded_endpoint, POLYGON_AREA),
    )
    if abs(cases[0]["time"] - cases[1]["time"]) > 1e-10:
        raise ValueError("Endpoints are not at identical native times")

    axis = np.linspace(*DOMAIN, PIXELS)
    xx, yy = np.meshgrid(axis, axis)
    query = np.column_stack((xx.ravel(), yy.ravel()))
    theta = np.arange(1024) * 2 * np.pi / 1024
    boundary = np.column_stack((3 + 2.75 * np.cos(theta),
                                3 + 2.75 * np.sin(theta)))
    inside = MplPath(boundary).contains_points(query).reshape(PIXELS, PIXELS)
    fig, panels = plt.subplots(2, 2, figsize=(10.8, 9.2),
                               constrained_layout=True)
    for col, case in enumerate(cases):
        nearest = cKDTree(case["xy"]).query(query, k=1, workers=1)[1]
        density = case["density"][nearest].reshape(PIXELS, PIXELS)
        relative = case["relative"][nearest].reshape(PIXELS, PIXELS)
        if col == 1:
            density = np.ma.masked_where(~inside, density)
            relative = np.where(inside, relative, np.nan)
        density_artist = panels[0, col].imshow(
            density, origin="lower", extent=(*DOMAIN, *DOMAIN),
            norm=LogNorm(vmin=.01, vmax=1.1, clip=True), cmap="viridis",
            interpolation="nearest")
        relative_artist = panels[1, col].imshow(
            np.ma.masked_invalid(relative), origin="lower",
            extent=(*DOMAIN, *DOMAIN),
            norm=SymLogNorm(linthresh=1e-4, linscale=1., vmin=-.1,
                            vmax=.1, base=10), cmap="coolwarm",
            interpolation="nearest")
        panels[0, col].set_title(f"{case['label']} ({case['count']:,} cells)")
        panels[1, col].set_title(f"Disk mean |ΔΣ/Σ₀| = {case['disk_abs']:.3g}")
    fig.colorbar(density_artist, ax=panels[0, :].tolist(), label="surface density Σ")
    fig.colorbar(relative_artist, ax=panels[1, :].tolist(), label="disk Σ/Σ₀ − 1")
    for panel in panels.flat:
        panel.set_aspect("equal")
        panel.set_xlabel("x")
        panel.set_ylabel("y")
        panel.set_xticks((1, 2, 3, 4, 5))
        panel.set_yticks((1, 2, 3, 4, 5))
    fig.suptitle(f"Native Kepler density at t={cases[0]['time']:.6g}; "
                 "common [-2,2] view and colors")
    fig.savefig(args.output, dpi=150)
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()

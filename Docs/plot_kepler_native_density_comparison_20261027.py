"""Plot same-time native W=0 Kepler density fields on one common scale."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
import numpy as np
from scipy.spatial import cKDTree

from plot_kepler_native_density_long_20261026 import complete_ordered


DOMAIN = (1., 5.)  # The paper's [-2,2] view, shifted to this run's center (3,3).
PIXELS = 720


def sample_nearest(xy, values, query):
    nearest = cKDTree(xy).query(query, k=1, workers=1)[1]
    return np.asarray(values)[nearest].reshape(PIXELS, PIXELS)


def load_case(label, initial, endpoint):
    p0, m0 = complete_ordered(initial)
    p1, m1 = complete_ordered(endpoint)
    if len(p0) != len(p1) or float(m0["t"]) != 0 or float(m1["t"]) <= 0:
        raise ValueError(f"{label}: mismatched initial/endpoint snapshots")
    if (float(m0["Lx"]), float(m0["Ly"]), float(m1["Lx"]), float(m1["Ly"])) != (6., 6., 6., 6.):
        raise ValueError(f"{label}: expected 6x6 Kepler box")
    if not np.array_equal(p0["mass"], p1["mass"]):
        raise ValueError(f"{label}: cell masses changed")
    if not (np.all(p0["w2"] == 0) and np.all(p1["w2"] == 0)):
        raise ValueError(f"{label}: nearest-site raster requires W=0")
    for p in (p0, p1):
        if not all(np.isfinite(p[key]).all() for key in ("x", "y", "mass", "volume", "den")):
            raise ValueError(f"{label}: nonfinite state")
        if not (np.all(p["mass"] > 0) and np.all(p["volume"] > 0) and np.all(p["den"] > 0)):
            raise ValueError(f"{label}: nonpositive state")
        if abs(float(np.sum(p["volume"], dtype=np.longdouble)) - 36.) > 1e-4:
            raise ValueError(f"{label}: cells do not tile the periodic box")
        if np.max(abs(p["den"].astype(float) /
                      (p["mass"].astype(float) / p["volume"].astype(float)) - 1.)) > 1e-6:
            raise ValueError(f"{label}: density differs from mass/volume")
    r0 = np.hypot(p0["x"] - 3., p0["y"] - 3.)
    disk = (r0 >= .5) & (r0 < 2.)
    mass = p0["mass"].astype(float)
    relative = p1["den"].astype(float) / p0["den"].astype(float) - 1.
    disk_error = float(np.sum(abs(relative[disk]) * mass[disk], dtype=np.longdouble) /
                       np.sum(mass[disk], dtype=np.longdouble))
    relative[~disk] = np.nan
    return {"label": label, "time": float(m1["t"]), "count": len(p1),
            "xy": np.column_stack((p1["x"], p1["y"])),
            "density": p1["den"].astype(float), "relative": relative,
            "disk_error": disk_error}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", nargs=3,
                        metavar=("LABEL", "INITIAL", "ENDPOINT"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--time-tolerance", type=float, default=1e-10,
                        help="maximum allowed endpoint-time difference (default: 1e-10)")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if not np.isfinite(args.time_tolerance) or not (0 <= args.time_tolerance <= 1e-6):
        raise ValueError("Time tolerance must be finite and between 0 and 1e-6")
    cases = [load_case(label, Path(initial), Path(endpoint))
             for label, initial, endpoint in args.case]
    if any(abs(case["time"] - cases[0]["time"]) > args.time_tolerance
           for case in cases[1:]):
        raise ValueError("All native endpoints must be at the same physical time")

    axis = np.linspace(*DOMAIN, PIXELS)
    xx, yy = np.meshgrid(axis, axis)
    query = np.column_stack((xx.ravel(), yy.ravel()))
    fig, panels = plt.subplots(2, len(cases), squeeze=False,
                               figsize=(4.8 * len(cases), 8.3), constrained_layout=True)
    for j, case in enumerate(cases):
        density = sample_nearest(case["xy"], case["density"], query)
        relative = sample_nearest(case["xy"], case["relative"], query)
        top = panels[0, j].imshow(density, origin="lower", extent=(*DOMAIN, *DOMAIN),
                                  norm=LogNorm(vmin=.01, vmax=1.1, clip=True),
                                  cmap="viridis", interpolation="nearest")
        bottom = panels[1, j].imshow(np.ma.masked_invalid(relative), origin="lower",
                                     extent=(*DOMAIN, *DOMAIN),
                                     norm=SymLogNorm(linthresh=1e-4, linscale=1.,
                                                     vmin=-.1, vmax=.1, base=10),
                                     cmap="coolwarm", interpolation="nearest")
        panels[0, j].set_title(f"{case['label']}  ({case['count']:,} cells)")
        panels[1, j].set_title(f"disk mean |ΔΣ/Σ₀| = {case['disk_error']:.3g}")
    fig.colorbar(top, ax=panels[0, :].tolist(), label="surface density Σ")
    fig.colorbar(bottom, ax=panels[1, :].tolist(), label="disk Σ/Σ₀ − 1")
    for panel in panels.flat:
        panel.set_aspect("equal")
        panel.set_xlabel("x")
        panel.set_ylabel("y")
        panel.set_xticks((1, 2, 3, 4, 5))
        panel.set_yticks((1, 2, 3, 4, 5))
    fig.suptitle(f"Native Voronoi density at t={cases[0]['time']:.6g}\n"
                 "GIZMO [-2,2] view; common color scales")
    fig.savefig(args.output, dpi=160)
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()

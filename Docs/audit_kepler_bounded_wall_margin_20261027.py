"""Read-only generator clearance from the exact static 1,024-gon wall."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


COUNT = 55351
SIDES = 1024
RADIUS = 2.75
STEP = 2 * np.pi / SIDES
APOTHEM = RADIUS * np.cos(STEP / 2)
THRESHOLDS = (.105, .11, .15, .2, .3)


def describe(path, initial_mass):
    time = snapshot_time(path, COUNT)
    state, _ = ordered(path, COUNT)
    x = state["x"].astype(float) - 3.
    y = state["y"].astype(float) - 3.
    radius = np.hypot(x, y)
    theta = np.mod(np.arctan2(y, x), 2 * np.pi)
    side = np.floor(theta / STEP).astype(np.int64) % SIDES
    angle = (side + .5) * STEP
    nx, ny = np.cos(angle), np.sin(angle)
    gap = APOTHEM - (x * nx + y * ny)
    vn = state["vx"].astype(float) * nx + state["vy"].astype(float) * ny
    vr = (x * state["vx"].astype(float) + y * state["vy"].astype(float)) / np.maximum(radius, 1e-30)
    closest = np.argsort(gap)[:10]
    nearest = int(closest[0])
    near = gap < .2
    outward = np.maximum(vn[near], 0.)
    rows = []
    for i in closest:
        i = int(i)
        rows.append({
            "id": i, "wall_side": int(side[i]),
            "radius": float(radius[i]), "gap": float(gap[i]),
            "wall_normal_speed": float(vn[i]),
            "radial_speed": float(vr[i]),
            "density": float(state["den"][i]),
            "pressure": float(state["pressure"][i]),
        })
    return {
        "snapshot": str(path.resolve()), "sha256": digest(path), "time": time,
        "fixed_mass_exact": bool(np.array_equal(state["mass"], initial_mass)),
        "W0_exact": bool(np.all(state["w2"] == 0)),
        "all_generators_inside_polygon": bool(np.all(gap > 0)),
        "minimum_wall_gap": float(gap[nearest]),
        "minimum_gap_site": rows[0],
        "ten_nearest_sites": rows,
        "sites_with_gap_below": {str(t): int(np.count_nonzero(gap < t)) for t in THRESHOLDS},
        "near_wall_gap_lt_0p2": {
            "sites": int(np.count_nonzero(near)),
            "maximum_outward_normal_speed": float(np.max(outward)) if len(outward) else 0.,
            "rms_wall_normal_speed": float(np.sqrt(np.mean(vn[near]**2))) if np.any(near) else 0.,
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("endpoints", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if abs(snapshot_time(args.initial, COUNT)) > 1e-12:
        raise ValueError("Initial state is not exact t=0")
    initial, _ = ordered(args.initial, COUNT)
    initial_mass = initial["mass"].copy()
    paths = (args.initial, *args.endpoints)
    series = [describe(path, initial_mass) for path in paths]
    if any(series[i]["time"] >= series[i+1]["time"] for i in range(len(series)-1)):
        raise ValueError("Snapshots must be in strictly increasing physical time")
    if not all(row["all_generators_inside_polygon"] and row["fixed_mass_exact"]
               and row["W0_exact"] for row in series):
        raise ValueError("Wall, fixed-mass, or W0 gate failed")
    result = {
        "initial": str(args.initial.resolve()),
        "polygon_sides": SIDES, "circumradius": RADIUS,
        "apothem": float(APOTHEM), "series": series,
        "scope": ("Generator clearance to the static polygon. Instantaneous "
                  "wall-normal velocity and gap/speed quotient do not predict "
                  "future impact or integrate wall torque."),
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

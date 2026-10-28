"""Track angular-momentum exchange between fixed initial-radius site groups."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


COUNT = 55351
EDGES = (0.0, 0.5, 0.75, 1.0, 2.0, 2.35, 2.6, float("inf"))


def momentum(state, mass):
    x = state["x"].astype(float) - 3.0
    y = state["y"].astype(float) - 3.0
    return mass * (x * state["vy"].astype(float) - y * state["vx"].astype(float))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("wall_series", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    wall = json.loads(args.wall_series.read_text())
    series = wall["series"]
    if len(series) < 2 or abs(series[0]["time"]) > 1e-12:
        raise ValueError("Expected a t=0 reference and later snapshots")
    reference, _ = ordered(Path(series[0]["snapshot"]), COUNT)
    mass = reference["mass"]
    radius = np.hypot(reference["x"].astype(float) - 3.0,
                      reference["y"].astype(float) - 3.0)
    masks = [((radius >= low) & (radius < high)) for low, high in zip(EDGES[:-1], EDGES[1:])]
    labels = [f"[{low:g},{high:g})" for low, high in zip(EDGES[:-1], EDGES[1:])]
    initial_l = momentum(reference, mass.astype(float))
    initial_total = float(np.sum(initial_l, dtype=np.longdouble))
    initial_by_band = [float(np.sum(initial_l[mask], dtype=np.longdouble)) for mask in masks]
    output = []
    for row in series:
        path = Path(row["snapshot"])
        if digest(path) != row["sha256"]:
            raise ValueError(f"Snapshot hash differs from wall audit: {path}")
        if abs(snapshot_time(path, COUNT) - row["time"]) > 1e-10:
            raise ValueError(f"Snapshot time differs from wall audit: {path}")
        state, _ = ordered(path, COUNT)
        if not np.array_equal(state["mass"], mass) or not np.all(state["w2"] == 0):
            raise ValueError(f"Fixed-mass W=0 contract failed: {path}")
        l = momentum(state, mass.astype(float))
        total = float(np.sum(l, dtype=np.longdouble))
        bands = []
        for label, mask, base in zip(labels, masks, initial_by_band):
            current = float(np.sum(l[mask], dtype=np.longdouble))
            bands.append({"initial_radius": label, "sites": int(mask.sum()),
                          "initial_angular_momentum": base,
                          "delta_angular_momentum": current - base,
                          "relative_change": current / base - 1 if base else None})
        output.append({"time": row["time"], "snapshot": str(path.resolve()),
                       "sha256": row["sha256"], "total_angular_momentum": total,
                       "total_relative_change": total / initial_total - 1,
                       "bands": bands})
    result = {"wall_series": str(args.wall_series.resolve()),
              "initial_total_angular_momentum": initial_total, "series": output,
              "scope": "Material groups selected by initial site radius. Group changes measure exchange, not a wall-torque ledger or a numerical error attribution."}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"t={output[-1]['time']:.6g} total relative L={output[-1]['total_relative_change']:+.9e}")
    for band in output[-1]["bands"]:
        print(band["initial_radius"], band["sites"],
              f"delta L={band['delta_angular_momentum']:+.9e}",
              f"relative={band['relative_change']:+.9e}")


if __name__ == "__main__":
    main()

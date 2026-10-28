"""Locate where the bounded and periodic q2 exteriors first diverge."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from audit_kepler_sparse_area_t0825_20261027 import ordered
from compare_kepler_bounded_periodic_disk_t0825_20261027 import (
    N_BOUNDED, N_PERIODIC, TIME, entropy_change, weighted_abs,
)
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


BINS = ((0.0, 0.5), (0.5, 2.0), (2.0, 2.35), (2.35, 2.5), (2.5, 2.65))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bounded_initial", type=Path)
    ap.add_argument("bounded_endpoint", type=Path)
    ap.add_argument("periodic_initial", type=Path)
    ap.add_argument("periodic_endpoint", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    paths = (args.bounded_initial, args.bounded_endpoint,
             args.periodic_initial, args.periodic_endpoint)
    for path, n, t in zip(paths, (N_BOUNDED, N_BOUNDED, N_PERIODIC, N_PERIODIC),
                          (0.0, TIME, 0.0, TIME)):
        if abs(snapshot_time(path, n) - t) > 1e-10:
            raise ValueError(f"Invalid native checkpoint: {path}")
    b0, _ = ordered(paths[0], N_BOUNDED)
    b1, _ = ordered(paths[1], N_BOUNDED)
    p0, _ = ordered(paths[2], N_PERIODIC)
    p1, _ = ordered(paths[3], N_PERIODIC)
    xy0 = np.column_stack((b0["x"].astype(float), b0["y"].astype(float)))
    xp0 = np.column_stack((p0["x"].astype(float), p0["y"].astype(float)))
    distance, match = cKDTree(xp0).query(xy0, k=1, workers=1)
    if np.max(distance) != 0 or len(np.unique(match)) != N_BOUNDED:
        raise ValueError("Retained initial generators do not match")
    radius = np.hypot(xy0[:, 0] - 3, xy0[:, 1] - 3)
    b_den = b1["den"].astype(float) / b0["den"].astype(float) - 1
    p_den = p1["den"][match].astype(float) / p0["den"][match].astype(float) - 1
    b_ent = entropy_change(b0, b1, np.arange(N_BOUNDED))
    p_ent = entropy_change(p0, p1, match)
    position_diff = np.hypot(b1["x"].astype(float) - p1["x"][match].astype(float),
                             b1["y"].astype(float) - p1["y"][match].astype(float))
    rows = []
    for lo, hi in BINS:
        use = (radius >= lo) & (radius < hi)
        mass = b0["mass"][use].astype(float)
        rows.append({
            "rlo": lo, "rhi": hi, "sites": int(np.sum(use)),
            "bounded_initial_mass": float(np.sum(mass, dtype=np.longdouble)),
            "initial_mass_relative_max": float(np.max(np.abs(
                b0["mass"][use].astype(float)
                / p0["mass"][match[use]].astype(float) - 1))),
            "bounded_density_change_abs_mass_weighted": weighted_abs(b_den[use], mass),
            "periodic_density_change_abs_mass_weighted": weighted_abs(p_den[use], mass),
            "density_change_difference_abs_mass_weighted": weighted_abs(
                b_den[use] - p_den[use], mass),
            "bounded_entropy_log_abs_mass_weighted": weighted_abs(b_ent[use], mass),
            "periodic_entropy_log_abs_mass_weighted": weighted_abs(p_ent[use], mass),
            "position_difference_abs_mass_weighted": weighted_abs(position_diff[use], mass),
            "position_difference_max": float(np.max(position_diff[use])),
        })
    result = {
        "time": TIME,
        "bounded_initial_sha256": digest(paths[0]),
        "bounded_endpoint_sha256": digest(paths[1]),
        "periodic_initial_sha256": digest(paths[2]),
        "periodic_endpoint_sha256": digest(paths[3]),
        "matched_retained_sites": N_BOUNDED,
        "radial_bins": rows,
        "scope": "Matched generators; initial outer cell masses and domains can differ; early-time boundary localization only",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

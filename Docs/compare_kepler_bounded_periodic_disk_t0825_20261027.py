"""Same-generator disk comparison of bounded and sparse periodic Kepler."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N_BOUNDED = 55351
N_PERIODIC = 55931
TIME = 0.825
GAMMA_MINUS_ONE = float(np.float32(1.666667)) - 1.0


def weighted_abs(x, mass):
    return float(np.sum(np.abs(x) * mass, dtype=np.longdouble)
                 / np.sum(mass, dtype=np.longdouble))


def entropy_change(initial, final, indices):
    return (np.log(final["ie"][indices].astype(float)
                   / initial["ie"][indices].astype(float))
            + GAMMA_MINUS_ONE
            * np.log(final["volume"][indices].astype(float)
                     / initial["volume"][indices].astype(float)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bounded_initial", type=Path)
    parser.add_argument("bounded_endpoint", type=Path)
    parser.add_argument("periodic_initial", type=Path)
    parser.add_argument("periodic_endpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    paths = (args.bounded_initial, args.bounded_endpoint,
             args.periodic_initial, args.periodic_endpoint)
    expected = (0.0, TIME, 0.0, TIME)
    counts = (N_BOUNDED, N_BOUNDED, N_PERIODIC, N_PERIODIC)
    for path, n, t in zip(paths, counts, expected):
        if abs(snapshot_time(path, n) - t) > 1e-10:
            raise ValueError(f"Incomplete or wrong-time snapshot: {path}")
    b0, _ = ordered(paths[0], N_BOUNDED)
    b1, _ = ordered(paths[1], N_BOUNDED)
    p0, _ = ordered(paths[2], N_PERIODIC)
    p1, _ = ordered(paths[3], N_PERIODIC)
    xy_b0 = np.column_stack((b0["x"].astype(float), b0["y"].astype(float)))
    xy_p0 = np.column_stack((p0["x"].astype(float), p0["y"].astype(float)))
    distance, match = cKDTree(xy_p0).query(xy_b0, k=1, workers=1)
    radius = np.hypot(xy_b0[:, 0] - 3, xy_b0[:, 1] - 3)
    disk = np.flatnonzero((radius >= 0.5) & (radius < 2.0))
    matched = match[disk]
    mass = b0["mass"][disk].astype(float)
    b_entropy = entropy_change(b0, b1, disk)
    p_entropy = entropy_change(p0, p1, matched)
    b_density = b1["den"][disk].astype(float) / b0["den"][disk].astype(float) - 1
    p_density = p1["den"][matched].astype(float) / p0["den"][matched].astype(float) - 1
    b_x = b1["x"][disk].astype(float)
    b_y = b1["y"][disk].astype(float)
    p_x = p1["x"][matched].astype(float)
    p_y = p1["y"][matched].astype(float)
    position_difference = np.hypot(b_x - p_x, b_y - p_y)
    velocity_difference = np.hypot(
        b1["vx"][disk].astype(float) - p1["vx"][matched].astype(float),
        b1["vy"][disk].astype(float) - p1["vy"][matched].astype(float))
    mass_initial_rel = np.abs(b0["mass"][disk].astype(float)
                              / p0["mass"][matched].astype(float) - 1)
    initial_area_rel = np.abs(b0["volume"][disk].astype(float)
                              / p0["volume"][matched].astype(float) - 1)
    checks = {
        "matched_original_generators": bool(np.max(distance[disk]) == 0),
        "unique_periodic_matches": bool(len(np.unique(matched)) == len(disk)),
        "initial_disk_mass_matches": bool(np.max(mass_initial_rel) < 1e-6),
        "initial_disk_areas_match": bool(np.max(initial_area_rel) < 1e-6),
        "fixed_mass_both": bool(np.array_equal(b0["mass"], b1["mass"])
                                and np.array_equal(p0["mass"], p1["mass"])),
        "W0_both": bool(np.all(b0["w2"] == 0) and np.all(b1["w2"] == 0)
                        and np.all(p0["w2"] == 0) and np.all(p1["w2"] == 0)),
    }
    result = {
        "time": TIME,
        "bounded_initial_sha256": digest(paths[0]),
        "bounded_endpoint_sha256": digest(paths[1]),
        "periodic_initial_sha256": digest(paths[2]),
        "periodic_endpoint_sha256": digest(paths[3]),
        "disk_cells": len(disk),
        "initial_generator_distance_max": float(np.max(distance[disk])),
        "initial_mass_relative_max": float(np.max(mass_initial_rel)),
        "initial_area_relative_max": float(np.max(initial_area_rel)),
        "bounded_density_change_abs_mass_weighted": weighted_abs(b_density, mass),
        "periodic_density_change_abs_mass_weighted": weighted_abs(p_density, mass),
        "matched_density_change_difference_abs_mass_weighted": weighted_abs(
            b_density - p_density, mass),
        "bounded_entropy_log_abs_mass_weighted": weighted_abs(b_entropy, mass),
        "periodic_entropy_log_abs_mass_weighted": weighted_abs(p_entropy, mass),
        "matched_entropy_change_difference_abs_mass_weighted": weighted_abs(
            b_entropy - p_entropy, mass),
        "matched_position_difference_abs_mass_weighted": weighted_abs(
            position_difference, mass),
        "matched_velocity_difference_abs_mass_weighted": weighted_abs(
            velocity_difference, mass),
        "checks": checks,
        "scope": "Same initial disk sites and native time; distinct exterior domain and wall boundary; early comparison only",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))
    if not all(checks.values()):
        raise SystemExit("Same-disk matching failed")


if __name__ == "__main__":
    main()

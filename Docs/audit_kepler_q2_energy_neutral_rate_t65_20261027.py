"""Screen a fixed-mass, pressure-work-neutral correction to Jacobi Wdot.

This only checks the single global power constraint and linearized volumes.
It does not establish a local physical weight law or long-time stability.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from export_kepler_q2_weight_rate_csr_t65_20261027 import finite_power_edges
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from laguerre_weight_rate_actual_t625_20261026 import face_volume_rate
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
DT = 0.00193983


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("accepted_receipt", type=Path)
    parser.add_argument("reference_rates", type=Path)
    parser.add_argument("case", choices=("aligned", "jitter"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads(args.accepted_receipt.read_text())
    if receipt["time"] != 6.5 or not all(receipt["checks"].values()):
        raise ValueError("Expected an accepted t=6.5 checkpoint")
    if digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]:
        raise ValueError("Checkpoint hash changed")
    p0, _ = ordered(Path(receipt["initial"]), N)
    p, _ = ordered(Path(receipt["endpoint"]), N)
    xy0 = np.column_stack((p0["x"].astype(float)-3, p0["y"].astype(float)-3))
    xy = np.column_stack((p["x"].astype(float)-3, p["y"].astype(float)-3))
    v = np.column_stack((p["vx"].astype(float), p["vy"].astype(float)))
    r0 = np.linalg.norm(xy0, axis=1)
    active = (r0 >= .5) & (r0 < 2.)
    inner = active & (r0 < .75)
    lookup = np.full(N, -1, dtype=int)
    lookup[active] = np.arange(active.sum())
    with np.load(args.reference_rates) as saved:
        if not np.array_equal(saved["active"].astype(bool), active):
            raise ValueError("Active set mismatch")
        rate = saved[args.case].copy()
    area, meta, triangles, centers = power_areas(
        xy, np.zeros(N), np.flatnonzero(active), lookup, return_geometry=True)
    if meta["hidden_active"] or np.any(area <= 0):
        raise ValueError("Invalid power geometry")
    geometric_rate, _, _ = face_volume_rate(
        xy, v, np.zeros(N), np.zeros(N), np.flatnonzero(active),
        lookup, triangles, centers)
    i, j, c1, c2 = finite_power_edges(triangles, centers, active)
    conductance = np.linalg.norm(c2-c1, axis=1) / (2*np.linalg.norm(xy[i]-xy[j], axis=1))
    pressure = p["pressure"].astype(float)
    gradient = np.zeros(N)
    pair_gradient = conductance*(pressure[i]-pressure[j])
    np.add.at(gradient, i, pair_gradient)
    np.add.at(gradient, j, -pair_gradient)
    gradient[~active] = 0
    original_power = float(np.dot(gradient, rate))
    projected = rate - gradient*(original_power/np.dot(gradient, gradient))
    projected[~active] = 0

    def measure(candidate):
        pair = conductance*(candidate[i]-candidate[j])
        dv = np.zeros(N)
        np.add.at(dv, i, pair)
        np.add.at(dv, j, -pair)
        residual = geometric_rate + dv[active]
        relative = np.abs(DT*residual/area)
        relative_inner = relative[inner[active]]
        next_weight = DT*candidate
        next_area, next_meta = power_areas(
            xy+DT*v, next_weight, np.flatnonzero(active), lookup)
        if next_meta["hidden_active"] or np.any(next_area <= 0):
            raise ValueError("Candidate hides a disk power cell")
        rebuilt = np.abs(next_area/area-1)
        face_guard = np.max(np.abs(next_weight[i]-next_weight[j])
                            / np.sum((xy[i]-xy[j])**2, axis=1))
        return {
            "pressure_weight_power": float(np.dot(gradient, candidate)),
            "area_change_p99_linearized": float(np.quantile(relative, .99)),
            "inner_area_change_p99_linearized": float(np.quantile(relative_inner, .99)),
            "area_change_max_linearized": float(relative.max()),
            "area_change_p99_rebuilt": float(np.quantile(rebuilt, .99)),
            "inner_area_change_p99_rebuilt": float(np.quantile(rebuilt[inner[active]], .99)),
            "area_change_max_rebuilt": float(rebuilt.max()),
            "max_weight_difference_over_pair_distance_squared": float(face_guard),
            "rate_abs_max": float(np.max(np.abs(candidate))),
        }

    result = {
        "case": args.case, "endpoint_sha256": receipt["endpoint_sha256"],
        "reference_rates_sha256": digest(args.reference_rates),
        "active_cells": int(active.sum()),
        "geometry_only": measure(np.zeros(N)),
        "jacobi16": measure(rate),
        "energy_neutral_projection": measure(projected),
        "relative_rate_change_l2": float(np.linalg.norm(projected-rate)/np.linalg.norm(rate)),
        "scope": "One accepted W0 t=6.5 state, linearized and independently rebuilt one-Euler-step power geometry. The algebraic global energy-neutral projection is only a candidate constraint; it does not prove local momentum/entropy conservation or a material-face model.",
    }
    if (abs(result["jacobi16"]["pressure_weight_power"]-original_power) > 1e-20
            or abs(result["energy_neutral_projection"]["pressure_weight_power"]) > 1e-18):
        raise ValueError("Energy-neutral projection check failed")
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

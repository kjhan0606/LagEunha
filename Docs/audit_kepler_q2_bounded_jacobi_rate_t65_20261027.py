"""Screen synchronization-free Jacobi weight-rate iterations against CG16."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_newton_volume_t6_20261026 import newton_jacobian
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from laguerre_weight_rate_actual_t625_20261026 import face_volume_rate
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
DT = .00193983
COUNTS = (8, 16, 32, 64)


def summary(rate, rhs, matrix, diagonal, target, inner):
    residual = rhs - matrix @ rate
    predicted = np.abs(DT*residual/target)
    return {"residual_l2_over_rhs": float(np.linalg.norm(residual)/np.linalg.norm(rhs)),
            "linearized_area_change_p99": float(np.quantile(predicted, .99)),
            "linearized_inner_p99": float(np.quantile(predicted[inner], .99)),
            "linearized_area_change_max": float(predicted.max()),
            "rate_abs_max": float(np.max(np.abs(rate))),
            "diagonal_min": float(diagonal.min())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accepted_receipt", type=Path)
    ap.add_argument("static_weights", type=Path)
    ap.add_argument("cg_rates", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    receipt = json.loads(args.accepted_receipt.read_text())
    if (receipt["time"] != 6.5 or len(receipt["checks"]) != 17
            or not all(receipt["checks"].values())):
        raise ValueError("Unaccepted t=6.5 native state")
    if digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]:
        raise ValueError("Accepted checkpoint changed")
    p0, _ = ordered(Path(receipt["initial"]), N)
    p, _ = ordered(Path(receipt["endpoint"]), N)
    xy0 = np.column_stack((p0["x"].astype(float)-3., p0["y"].astype(float)-3.))
    xy = np.column_stack((p["x"].astype(float)-3., p["y"].astype(float)-3.))
    velocity = np.column_stack((p["vx"].astype(float), p["vy"].astype(float)))
    r0 = np.linalg.norm(xy0, axis=1)
    active = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    lookup = np.full(N, -1, dtype=int)
    lookup[active] = np.arange(len(active))
    with np.load(args.static_weights) as saved:
        if not np.array_equal(saved["active"], active):
            raise ValueError("Static active set changed")
        weight = np.zeros(N)
        weight[active] = saved["weight"]
    with np.load(args.cg_rates) as saved:
        if not np.array_equal(saved["active"], active):
            raise ValueError("CG active set changed")
        cg16 = saved["rate_16"].copy()
    target = p0["volume"][active].astype(float)
    area, meta, triangles, centers = power_areas(
        xy, weight, active, lookup, return_geometry=True)
    if meta["hidden_active"] or np.max(np.abs(area/target-1)) > 1e-7:
        raise ValueError("Static power volume mismatch")
    jac, diagonal = newton_jacobian(xy, lookup, len(active), triangles, centers)
    fixed_rate, _, _ = face_volume_rate(
        xy, velocity, weight, np.zeros(N), active, lookup, triangles, centers)
    rhs = -fixed_rate
    inner = (r0[active] >= .5) & (r0[active] < .75)
    trials = {}
    best = None
    for omega in (.5, 2./3., .8, 1.):
        rate = np.zeros(len(active))
        row = {}
        for iteration in range(1, max(COUNTS)+1):
            rate += omega*(rhs-jac@rate)/diagonal
            if iteration in COUNTS:
                item = summary(rate, rhs, jac, diagonal, target, inner)
                row[str(iteration)] = item
                if iteration == 16 and (best is None or
                        item["linearized_inner_p99"] < best[0]):
                    best = (item["linearized_inner_p99"], omega, rate.copy())
        trials[str(omega)] = row
    cg_summary = summary(cg16, rhs, jac, diagonal, target, inner)
    best_rate = best[2]
    next_xy = xy+DT*velocity
    predicted_weight = weight.copy()
    predicted_weight[active] += DT*best_rate
    next_area, next_meta = power_areas(next_xy, predicted_weight, active, lookup)
    if next_meta["hidden_active"] or np.any(next_area <= 0):
        raise ValueError("Best Jacobi16 hides a power cell")
    actual = np.abs(next_area/target-1)
    result = {
        "endpoint_sha256": receipt["endpoint_sha256"],
        "static_weights_sha256": digest(args.static_weights),
        "cg_rates_sha256": digest(args.cg_rates),
        "active_cells": len(active),
        "target": "J Wdot = -fixed-W power-cell area rate",
        "step": DT,
        "cg16": cg_summary,
        "jacobi": trials,
        "best_jacobi16_omega_by_inner_p99": best[1],
        "best_jacobi16_rebuilt_area": {
            "whole_p99": float(np.quantile(actual, .99)),
            "inner_p99": float(np.quantile(actual[inner], .99)),
            "max": float(actual.max()),
        },
        "scope": "One accepted t=6.5 geometry and one Euler displacement, no native hydro or GPU timing. Jacobi screening against an existing 16-step PCG reference; the best relaxation was chosen using this same state and needs independent validation.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"best_omega": best[1],
                      "cg16_inner_p99_linear": cg_summary["linearized_inner_p99"],
                      "jacobi16_inner_p99_linear": best[0],
                      "jacobi16_inner_p99_rebuilt": result["best_jacobi16_rebuilt_area"]["inner_p99"]}, indent=2))


if __name__ == "__main__":
    main()

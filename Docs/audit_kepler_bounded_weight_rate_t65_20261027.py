"""One-stage geometric Laguerre weight-rate solve at native bounded t=6.5."""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator, cg

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_bounded_weight_feasibility_t65_20261027 import metrics
from laguerre_newton_volume_t6_20261026 import newton_jacobian
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from laguerre_weight_rate_actual_t625_20261026 import face_volume_rate
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
STEP = 0.00193983
COUNTS = (8, 16, 32, 64, 128, 256, 512, 1024)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accepted_receipt", type=Path)
    ap.add_argument("static_weights", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--rates", type=Path, required=True)
    args = ap.parse_args()
    started = time.monotonic()
    receipt = json.loads(args.accepted_receipt.read_text())
    if (receipt["time"] != 6.5 or len(receipt["checks"]) != 17
            or not all(receipt["checks"].values())):
        raise ValueError("Expected accepted native t=6.5 receipt")
    p0, _ = ordered(Path(receipt["initial"]), N)
    p, _ = ordered(Path(receipt["endpoint"]), N)
    if (digest(Path(receipt["initial"])) != receipt["initial_sha256"]
            or digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]):
        raise ValueError("Native snapshot hash differs from receipt")
    xy0 = np.column_stack((p0["x"].astype(float)-3., p0["y"].astype(float)-3.))
    xy = np.column_stack((p["x"].astype(float)-3., p["y"].astype(float)-3.))
    velocity = np.column_stack((p["vx"].astype(float), p["vy"].astype(float)))
    r0 = np.linalg.norm(xy0, axis=1)
    active = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    lookup = np.full(N, -1, dtype=int)
    lookup[active] = np.arange(len(active))
    with np.load(args.static_weights) as saved:
        if not np.array_equal(saved["active"], active):
            raise ValueError("Static weighted active set differs")
        local_weight = saved["weight"].copy()
    weight = np.zeros(N)
    weight[active] = local_weight
    target = p0["volume"][active].astype(float)
    area, geom, triangles, centers = power_areas(
        xy, weight, active, lookup, return_geometry=True)
    if geom["hidden_active"] or np.any(area <= 0) or np.max(np.abs(area/target-1)) > 1e-7:
        raise ValueError("Static volume solution was not restored")
    jac, diagonal = newton_jacobian(xy, lookup, len(active), triangles, centers)
    fixed_rate, _, nfaces = face_volume_rate(
        xy, velocity, weight, np.zeros(N), active, lookup, triangles, centers)
    h = 1e-5
    plus, plus_meta = power_areas(xy+h*velocity, weight, active, lookup)
    minus, minus_meta = power_areas(xy-h*velocity, weight, active, lookup)
    if plus_meta["hidden_active"] or minus_meta["hidden_active"]:
        raise ValueError("Area-rate finite difference hid an active cell")
    finite_difference = (plus-minus)/(2*h)
    face_fd_relative_l2 = float(np.linalg.norm(fixed_rate-finite_difference)
                                / np.linalg.norm(finite_difference))
    rhs = -fixed_rate
    snapshots = {}
    iterations = [0]

    def capture(vector):
        iterations[0] += 1
        if iterations[0] in COUNTS:
            snapshots[iterations[0]] = vector.copy()

    rate, info = cg(jac, rhs,
                    M=LinearOperator(jac.shape, matvec=lambda v: v/diagonal),
                    rtol=1e-10, atol=0., maxiter=2000, callback=capture)
    if info != 0:
        raise RuntimeError(f"Weight-rate CG failed after {iterations[0]} iterations: {info}")
    all_rate = np.zeros(N)
    all_rate[active] = rate
    _, weighted_face_rate, _ = face_volume_rate(
        xy, velocity, weight, all_rate, active, lookup, triangles, centers)
    next_xy = xy + STEP*velocity
    fixed_area, fixed_meta = power_areas(next_xy, weight, active, lookup)
    if fixed_meta["hidden_active"] or np.any(fixed_area <= 0):
        raise ValueError("Fixed-weight predicted diagram invalid")
    trials = {}
    trial_vectors = {**snapshots, iterations[0]: rate}
    for count in sorted(trial_vectors):
        local_rate = trial_vectors[count]
        predicted = weight.copy()
        predicted[active] += STEP*local_rate
        next_area, next_meta = power_areas(next_xy, predicted, active, lookup)
        if next_meta["hidden_active"] or np.any(~np.isfinite(next_area)) or np.any(next_area <= 0):
            trials[str(count)] = {"valid_power_cells": False}
            continue
        trials[str(count)] = {"valid_power_cells": True,
                              "cg_relative_residual": float(
                                  np.linalg.norm(rhs-jac@local_rate)/np.linalg.norm(rhs)),
                              "area_change_from_current": metrics(next_area, area, r0[active]),
                              "area_error_to_initial_target": metrics(next_area, target, r0[active])}
    result = {"accepted_receipt": str(args.accepted_receipt.resolve()),
              "endpoint_sha256": receipt["endpoint_sha256"],
              "static_weights_sha256": digest(args.static_weights),
              "time": 6.5, "active_disk_sites": len(active), "face_pairs_touching_active": nfaces,
              "geometry": geom, "step": STEP,
              "fixed_face_rate_vs_finite_difference_relative_l2": face_fd_relative_l2,
              "cold_cg_iterations": iterations[0],
              "full_weighted_face_rate_over_fixed_rate_l2": float(
                  np.linalg.norm(weighted_face_rate)/np.linalg.norm(fixed_rate)),
              "weight_rate_abs_quantiles": np.quantile(np.abs(rate), [.5,.95,.99,1.]).tolist(),
              "fixed_weight_area_change": metrics(fixed_area, area, r0[active]),
              "truncated_rate_trials": trials,
              "wall_seconds": time.monotonic()-started,
              "scope": "One forward-Euler geometric prediction from the actual t=6.5 state. No native weighted face force/work, hydro integration, GPU timing or long-time stability proof."}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    with args.rates.open("xb") as stream:
        np.savez_compressed(stream, active=active, weight_rate=rate,
                            **{f"rate_{key}": value for key, value in snapshots.items()})
    print("cg", iterations[0], "fd_l2", face_fd_relative_l2,
          "fixed_inner_p99", result["fixed_weight_area_change"]["inner_p99_relative"],
          "full_inner_p99", trials[str(iterations[0])]["area_change_from_current"]["inner_p99_relative"],
          "seconds", result["wall_seconds"])


if __name__ == "__main__":
    main()

"""Static Laguerre target-volume feasibility at an accepted bounded t=6.5 state.

This is a geometry reference, not a hydrodynamic weight evolution scheme.
"""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator, cg

from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_newton_volume_t6_20261026 import newton_jacobian
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N = 55351


def metrics(areas, target, initial_radius):
    delta = np.abs(areas / target - 1.)
    inner = (initial_radius >= .5) & (initial_radius < .75)
    return {"median_relative": float(np.median(delta)),
            "p95_relative": float(np.quantile(delta, .95)),
            "p99_relative": float(np.quantile(delta, .99)),
            "max_relative": float(np.max(delta)),
            "inner_median_relative": float(np.median(delta[inner])),
            "inner_p99_relative": float(np.quantile(delta[inner], .99)),
            "inner_max_relative": float(np.max(delta[inner]))}


def guard(xy, weights, triangles, lookup):
    edges = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]],
                       triangles[:, [2, 0]]))
    edges.sort(axis=1)
    edges = np.unique(edges, axis=0)
    touching = (lookup[edges[:, 0]] >= 0) | (lookup[edges[:, 1]] >= 0)
    edges = edges[touching]
    i, j = edges[:, 0], edges[:, 1]
    d2 = np.sum((xy[i]-xy[j])**2, axis=1)
    ratio = np.abs(weights[i]-weights[j]) / d2
    worst = np.argsort(ratio)[-10:][::-1]
    return {"faces_touching_active": len(edges),
            "max_absolute_weight_difference_over_d2": float(np.max(ratio)),
            "faces_over_0p8": int(np.count_nonzero(ratio > .8)),
            "faces_over_1": int(np.count_nonzero(ratio > 1.)),
            "worst_faces": [{"i": int(i[k]), "j": int(j[k]),
                             "ratio": float(ratio[k])} for k in worst]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("accepted_receipt", type=Path)
    parser.add_argument("--time", type=float, default=6.5)
    parser.add_argument("--expected-failed-check", action="append", default=[],
                        help="Exploratory geometry only: require exactly this failed native gate")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--weights", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    receipt = json.loads(args.accepted_receipt.read_text())
    failed = {key for key, passed in receipt["checks"].items() if not passed}
    if (receipt["time"] != args.time or len(receipt["checks"]) != 17
            or failed != set(args.expected_failed_check)):
        raise ValueError(f"Unexpected native time or failed checks: {sorted(failed)}")
    initial, endpoint = Path(receipt["initial"]), Path(receipt["endpoint"])
    if (digest(initial) != receipt["initial_sha256"]
            or digest(endpoint) != receipt["endpoint_sha256"]
            or abs(snapshot_time(initial, N)) > 1e-12
            or abs(snapshot_time(endpoint, N) - args.time) > 1e-10):
        raise ValueError("Native snapshots changed")
    p0, _ = ordered(initial, N)
    p, _ = ordered(endpoint, N)
    xy0 = np.column_stack((p0["x"].astype(float)-3., p0["y"].astype(float)-3.))
    xy = np.column_stack((p["x"].astype(float)-3., p["y"].astype(float)-3.))
    r0 = np.linalg.norm(xy0, axis=1)
    active = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    if len(active) != 49038:
        raise ValueError("Wrong active disk population")
    lookup = np.full(N, -1, dtype=int)
    lookup[active] = np.arange(len(active))
    target = p0["volume"][active].astype(float)
    weights = np.zeros(N)
    history = []
    initial_area_error = None
    converged = False
    for iteration in range(8):
        areas, geom, triangles, centers = power_areas(
            xy, weights, active, lookup, return_geometry=True)
        if geom["hidden_active"] or np.any(~np.isfinite(areas)) or np.any(areas <= 0):
            raise ValueError(f"Invalid active power cell at iteration {iteration}")
        row = {"iteration": iteration, "geometry": geom,
               "target_error": metrics(areas, target, r0[active])}
        history.append(row)
        if iteration == 0:
            initial_area_error = np.abs(areas / p["volume"][active].astype(float) - 1.)
            if np.max(initial_area_error) > 1e-5:
                raise ValueError("Reference power area differs from native Voronoi area")
        if row["target_error"]["max_relative"] < 1e-7:
            converged = True
            break
        jac, diagonal = newton_jacobian(xy, lookup, len(active), triangles, centers)
        if not np.all(diagonal > 0):
            raise ValueError("Nonpositive volume Jacobian diagonal")
        iterations = [0]
        correction, info = cg(
            jac, target-areas,
            M=LinearOperator(jac.shape, matvec=lambda v: v/diagonal),
            rtol=1e-9, atol=0, maxiter=2000,
            callback=lambda _: iterations.__setitem__(0, iterations[0]+1))
        row["cg_iterations"] = iterations[0]
        row["cg_info"] = int(info)
        if info != 0:
            break
        accepted = False
        current_norm = float(np.linalg.norm((areas-target)/target))
        for factor in (1., .5, .25, .125, .0625, .03125):
            candidate = weights.copy()
            candidate[active] += factor*correction
            next_area, next_geom = power_areas(xy, candidate, active, lookup)
            if next_geom["hidden_active"] or np.any(~np.isfinite(next_area)) or np.any(next_area <= 0):
                continue
            next_norm = float(np.linalg.norm((next_area-target)/target))
            if next_norm < current_norm:
                weights = candidate
                row["accepted_factor"] = factor
                accepted = True
                break
        if not accepted:
            row["accepted_factor"] = None
            break
    final_area, final_geom, final_triangles, _ = power_areas(
        xy, weights, active, lookup, return_geometry=True)
    if final_geom["hidden_active"] or np.any(final_area <= 0):
        raise ValueError("Invalid final power diagram")
    final = metrics(final_area, target, r0[active])
    result = {"accepted_receipt": str(args.accepted_receipt.resolve()),
              "initial_sha256": receipt["initial_sha256"],
              "endpoint_sha256": receipt["endpoint_sha256"],
              "time": args.time, "native_failed_checks": sorted(failed),
              "sites": N, "active_disk_sites": len(active),
              "target": "Each disk site's own native t=0 Voronoi volume; core and exterior weights fixed at zero",
              "native_voronoi_area_relative_max": float(np.max(initial_area_error)),
              "history": history, "converged_to_1e_minus_7": converged,
              "final_target_error": final,
              "weight_abs_quantiles": np.quantile(np.abs(weights[active]), [.5,.95,.99,1.]).tolist(),
              "weight_min": float(np.min(weights[active])),
              "weight_max": float(np.max(weights[active])),
              "native_face_guard": guard(xy, weights, final_triangles, lookup),
              "wall_seconds": time.monotonic()-started,
              "scope": "Static actual-position regular-triangulation geometry only. Listed native failed checks remain failures. No native GPU weighted tessellation, weight-rate integration, force/work/mass conservation or long-time result."}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    with args.weights.open("xb") as stream:
        np.savez_compressed(stream, active=active, weight=weights[active])
    print("converged", converged, "max_relative", final["max_relative"],
          "guard_over_0p8", result["native_face_guard"]["faces_over_0p8"],
          "wall_seconds", result["wall_seconds"])


if __name__ == "__main__":
    main()

"""Compare local face-velocity fitting with volume-constrained Wdot.

Read-only geometry/linear-solve diagnostic on accepted bounded q2 t=6.5 states.
The fit target is the native pair's mean velocity, with no analytic orbit input.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import LinearOperator, cg

from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
DT = 0.00193983
EPS = 0.01175


def touching_faces(triangles, centers, active_flag):
    edges = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]],
                       triangles[:, [2, 0]]))
    tid = np.tile(np.arange(len(triangles)), 3)
    edges.sort(axis=1)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    edges, tid = edges[order], tid[order]
    starts = np.r_[0, 1 + np.flatnonzero(np.any(edges[1:] != edges[:-1], axis=1))]
    counts = np.diff(np.r_[starts, len(edges)])
    touching = active_flag[edges[starts, 0]] | active_flag[edges[starts, 1]]
    if np.any(counts[touching] != 2):
        raise ValueError("A face touching the disk is not finite and paired")
    starts = starts[touching]
    i, j = edges[starts, 0], edges[starts, 1]
    c1, c2 = centers[tid[starts]], centers[tid[starts + 1]]
    length = np.linalg.norm(c2 - c1, axis=1)
    if not np.all(np.isfinite(length) & (length > 0)):
        raise ValueError("Invalid face length")
    return i, j, .5*(c1+c2), length


def face_metrics(slip, length, disk_disk):
    def one(mask):
        l, s = length[mask], slip[mask]
        return {"faces": int(len(s)),
                "length_weighted_mean_abs": float(np.dot(l, np.abs(s))/l.sum()),
                "length_weighted_rms": float(np.sqrt(np.dot(l, s*s)/l.sum()))}
    return {"all_touching_disk": one(np.ones(len(slip), dtype=bool)),
            "disk_disk": one(disk_disk)}


def volume_metrics(rate, area, inner):
    relative = np.abs(DT*rate/area)
    return {"median": float(np.median(relative)),
            "p99": float(np.quantile(relative, .99)),
            "inner_p99": float(np.quantile(relative[inner], .99)),
            "max": float(np.max(relative))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accepted_receipt", type=Path)
    ap.add_argument("static_weights", type=Path)
    ap.add_argument("volume_rates", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    receipt = json.loads(args.accepted_receipt.read_text())
    if (receipt["time"] != 6.5 or len(receipt["checks"]) != 17
            or not all(receipt["checks"].values())):
        raise ValueError("Unaccepted native t=6.5 state")
    if digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]:
        raise ValueError("Native checkpoint hash changed")
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
            raise ValueError("Static weight active set changed")
        weight = np.zeros(N)
        weight[active] = saved["weight"]
    with np.load(args.volume_rates) as saved:
        if not np.array_equal(saved["active"], active):
            raise ValueError("Volume-rate active set changed")
        volume16 = np.zeros(N)
        volume16[active] = saved["rate_16"]
    area, meta, triangles, centers = power_areas(
        xy, weight, active, lookup, return_geometry=True)
    if meta["hidden_active"] or np.any(~np.isfinite(area) | (area <= 0)):
        raise ValueError("Invalid weighted disk geometry")
    i, j, c, length = touching_faces(triangles, centers, lookup >= 0)
    dr = xy[j]-xy[i]
    distance = np.linalg.norm(dr, axis=1)
    normal = dr/distance[:, None]
    fixed = (np.einsum("ij,ij->i", xy[j], velocity[j])
             - np.einsum("ij,ij->i", xy[i], velocity[i])
             - np.einsum("ij,ij->i", velocity[j]-velocity[i], c))/distance
    local = np.einsum("ij,ij->i", .5*(velocity[i]+velocity[j]), normal)
    material = np.column_stack((-c[:, 1], c[:, 0]))
    material *= (np.einsum("ij,ij->i", c, c)+EPS*EPS)[:, None]**(-.75)
    analytic = np.einsum("ij,ij->i", material, normal)
    disk_disk = (lookup[i] >= 0) & (lookup[j] >= 0)

    face_id = np.arange(len(i))
    valid_i = lookup[i] >= 0
    valid_j = lookup[j] >= 0
    incidence = coo_matrix((np.r_[np.ones(valid_i.sum()), -np.ones(valid_j.sum())],
                            (np.r_[face_id[valid_i], face_id[valid_j]],
                             np.r_[lookup[i[valid_i]], lookup[j[valid_j]]])),
                           shape=(len(i), len(active))).tocsr()
    conductance = length/(4*distance*distance)
    matrix = (incidence.T @ incidence.multiply(conductance[:, None])).tocsr()
    rhs = incidence.T @ (length*(local-fixed)/(2*distance))
    diagonal = matrix.diagonal()
    if np.any(diagonal <= 0):
        raise ValueError("Singular active-cell diagonal")
    preconditioner = LinearOperator(matrix.shape, matvec=lambda x: x/diagonal)
    snapshots = {}
    iteration = [0]
    def record(x):
        iteration[0] += 1
        if iteration[0] in (8, 16, 32):
            snapshots[iteration[0]] = x.copy()
    solved, info = cg(matrix, rhs, M=preconditioner, rtol=1e-14,
                      atol=0., maxiter=32, callback=record)
    if info < 0:
        raise RuntimeError(f"Material-face CG broke down: {info}")
    for count in (8, 16, 32):
        if count not in snapshots:
            snapshots[count] = solved.copy()

    inner = (r0[active] >= .5) & (r0[active] < .75)
    all_cases = {"frozen": np.zeros(len(active)),
                 "volume16": volume16[active],
                 **{f"material{count}": snapshots[count] for count in (8, 16, 32)}}
    comparisons = {}
    for name, rate in all_cases.items():
        moving = fixed + incidence @ rate/(2*distance)
        vdot = incidence.T @ (length*moving)
        rate_global = np.zeros(N)
        rate_global[active] = rate
        guard = np.abs((weight[i]-weight[j]
                        +DT*(rate_global[i]-rate_global[j]))/(distance*distance))
        comparisons[name] = {
            "slip_to_local_velocity": face_metrics(moving-local, length, disk_disk),
            "slip_to_analytic_kepler": face_metrics(moving-analytic, length, disk_disk),
            "linearized_area_change": volume_metrics(vdot, area, inner),
            "face_guard_max_ratio": float(guard.max()),
            "faces_guard_over_0p8": int(np.count_nonzero(guard > .8)),
        }
    result = {
        "accepted_receipt": str(args.accepted_receipt.resolve()),
        "endpoint_sha256": receipt["endpoint_sha256"],
        "static_weights_sha256": digest(args.static_weights),
        "volume_rates_sha256": digest(args.volume_rates),
        "active_disk_cells": len(active),
        "touching_faces": len(i),
        "disk_disk_faces": int(disk_disk.sum()),
        "target": "Mean of the two native generator velocities, projected on each face normal; no analytic orbit in the solve",
        "material_fit_cg_iterations": iteration[0],
        "material_fit_cg_info": int(info),
        "material_fit_initial_rhs_l2": float(np.linalg.norm(rhs)),
        "material_fit_32_residual_over_rhs": float(np.linalg.norm(matrix@solved-rhs)/np.linalg.norm(rhs)),
        "cases": comparisons,
        "scope": "Read-only t=6.5 linearized face-rate comparison; no power diagram rebuilt at t+dt and no native weighted hydro evolution. The local gas velocity is a first-order pair proxy, not a reconstructed Riemann contact state.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"faces": len(i), "cg": iteration[0],
                      "cases": {name: {"local_slip_rms": d["slip_to_local_velocity"]["disk_disk"]["length_weighted_rms"],
                                       "analytic_slip_rms": d["slip_to_analytic_kepler"]["disk_disk"]["length_weighted_rms"],
                                       "area_p99": d["linearized_area_change"]["p99"]}
                                for name, d in comparisons.items()}}, indent=2))


if __name__ == "__main__":
    main()

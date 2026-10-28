"""Evaluate the missing pressure work of a prescribed Laguerre weight rate.

This is a frozen-state algebraic diagnostic, not a weighted hydro evolution.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_newton_volume_t6_20261026 import newton_jacobian
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
DT = 0.00193983
GAMMA = float(np.float32(1.666667))


def faces_touching_active(triangles, centers, active_flag):
    edges = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]],
                       triangles[:, [2, 0]]))
    tids = np.tile(np.arange(len(triangles)), 3)
    edges.sort(axis=1)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    edges, tids = edges[order], tids[order]
    starts = np.r_[0, 1+np.flatnonzero(np.any(edges[1:] != edges[:-1], axis=1))]
    counts = np.diff(np.r_[starts, len(edges)])
    touch = active_flag[edges[starts, 0]] | active_flag[edges[starts, 1]]
    if np.any(counts[touch] != 2):
        raise ValueError("Non-finite face touches active disk")
    starts = starts[touch]
    i, j = edges[starts, 0], edges[starts, 1]
    length = np.linalg.norm(centers[tids[starts]]-centers[tids[starts+1]], axis=1)
    return i, j, length


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
        raise ValueError("Native t=6.5 receipt is not accepted")
    for key, hash_key in (("initial", "initial_sha256"),
                          ("endpoint", "endpoint_sha256")):
        if digest(Path(receipt[key])) != receipt[hash_key]:
            raise ValueError(f"Changed native {key} snapshot")
    p0, _ = ordered(Path(receipt["initial"]), N)
    p, _ = ordered(Path(receipt["endpoint"]), N)
    xy0 = np.column_stack((p0["x"].astype(float)-3., p0["y"].astype(float)-3.))
    xy = np.column_stack((p["x"].astype(float)-3., p["y"].astype(float)-3.))
    r0 = np.linalg.norm(xy0, axis=1)
    active = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    lookup = np.full(N, -1, dtype=int)
    lookup[active] = np.arange(len(active))
    with np.load(args.static_weights) as saved:
        if not np.array_equal(saved["active"], active):
            raise ValueError("Static-weight active IDs differ")
        local_weight = saved["weight"].copy()
    with np.load(args.volume_rates) as saved:
        if not np.array_equal(saved["active"], active):
            raise ValueError("Weight-rate active IDs differ")
        local_rate = saved["weight_rate"].copy()
    weight, rate = np.zeros(N), np.zeros(N)
    weight[active], rate[active] = local_weight, local_rate
    area, meta, triangles, centers = power_areas(
        xy, weight, active, lookup, return_geometry=True)
    if (meta["hidden_active"] or np.any(area <= 0)
            or np.max(np.abs(area/p0["volume"][active]-1)) > 1e-7):
        raise ValueError("Static target-volume geometry changed")
    i, j, length = faces_touching_active(triangles, centers, lookup >= 0)
    distance = np.linalg.norm(xy[j]-xy[i], axis=1)
    conductance = length/(2*distance)
    delta_rate = rate[i]-rate[j]
    rate_volume_i = np.zeros(N)
    np.add.at(rate_volume_i, i, conductance*delta_rate)
    np.add.at(rate_volume_i, j, -conductance*delta_rate)
    jac, _ = newton_jacobian(xy, lookup, len(active), triangles, centers)
    jac_error = np.linalg.norm(rate_volume_i[active]-jac@local_rate)
    jac_scale = np.linalg.norm(rate_volume_i[active])
    if jac_error > 1e-10*jac_scale:
        raise ValueError("Face and Jacobian weight-volume rates disagree")
    pressure_native = p["pressure"].astype(float)
    if np.any(~np.isfinite(pressure_native)) or np.any(pressure_native <= 0):
        raise ValueError("Native pressure is invalid")
    pressure_reweighted = pressure_native.copy()
    pressure_reweighted[active] *= (
        p["volume"][active].astype(float)/area)**GAMMA
    mass = p["mass"].astype(float)
    kinetic = .5*np.sum(mass*(p["vx"].astype(float)**2+
                                p["vy"].astype(float)**2), dtype=np.longdouble)
    internal = np.sum(p["ie"].astype(float), dtype=np.longdouble)
    result = {"accepted_receipt": str(args.accepted_receipt.resolve()),
              "endpoint_sha256": receipt["endpoint_sha256"],
              "static_weights_sha256": digest(args.static_weights),
              "volume_rates_sha256": digest(args.volume_rates),
              "time": 6.5, "one_step_dt": DT,
              "disk_touching_faces": len(i),
              "disk_disk_faces": int(np.count_nonzero((lookup[i] >= 0) & (lookup[j] >= 0))),
              "jacobian_face_relative_l2_error": float(jac_error/jac_scale),
              "sum_weight_induced_volume_rate_all_cells": float(np.sum(rate_volume_i, dtype=np.longdouble)),
              "reference_kinetic_energy": float(kinetic),
              "reference_internal_energy": float(internal),
              "pressure_model": "Native t=6.5 pressure; optional disk-cell adiabatic correction to solved static power area. Inactive-cell pressure remains native, so corrected model is approximate at the disk boundary.",
              "scope": "Frozen-state pressure-weight work identity only; no native weighted force, flux, RK4, or long-time result."}
    for label, pressure in (("native_pressure", pressure_native),
                            ("reweighted_disk_pressure_proxy", pressure_reweighted)):
        face_work = conductance*(pressure[i]-pressure[j])*delta_rate
        cell_work = np.sum(pressure*rate_volume_i, dtype=np.longdouble)
        face_sum = np.sum(face_work, dtype=np.longdouble)
        face_absolute_sum = np.sum(np.abs(face_work), dtype=np.longdouble)
        if abs(float(cell_work-face_sum)) > 1e-12*max(float(face_absolute_sum),1e-30):
            raise ValueError("Face/cell pressure-weight work identity failed")
        result[label] = {
            "pressure_min": float(np.min(pressure)),
            "pressure_max": float(np.max(pressure)),
            "signed_pressure_weight_work_rate": float(face_sum),
            "absolute_face_work_rate_sum": float(face_absolute_sum),
            "one_step_signed_work_over_kinetic": float(DT*face_sum/kinetic),
            "one_step_absolute_face_work_over_kinetic": float(DT*face_absolute_sum/kinetic),
            "face_cell_identity_absolute_error": float(abs(cell_work-face_sum))}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: {"signed": result[k]["signed_pressure_weight_work_rate"],
                          "absolute": result[k]["absolute_face_work_rate_sum"]}
                      for k in ("native_pressure", "reweighted_disk_pressure_proxy")}))


if __name__ == "__main__":
    main()

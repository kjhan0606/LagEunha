"""Least-squares Laguerre weight rate for disk-face material motion at t=6.5."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import LinearOperator, cg

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_bounded_weight_feasibility_t65_20261027 import metrics as area_metrics
from audit_kepler_bounded_weighted_face_slip_t65_20261027 import finite_faces, metrics as slip_metrics
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from laguerre_weight_rate_actual_t625_20261026 import face_volume_rate
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
EPS = 0.01175
STEP = 0.00193983


def grouped(length, slip, same):
    return {"all_disk_disk": slip_metrics(length, slip),
            "same_initial_shell": slip_metrics(length[same], slip[same]),
            "different_initial_shell": slip_metrics(length[~same], slip[~same])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accepted_receipt", type=Path)
    ap.add_argument("static_weights", type=Path)
    ap.add_argument("volume_rates", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--rates", type=Path, required=True)
    args = ap.parse_args()
    receipt = json.loads(args.accepted_receipt.read_text())
    if (receipt["time"] != 6.5 or len(receipt["checks"]) != 17
            or not all(receipt["checks"].values())
            or digest(Path(receipt["initial"])) != receipt["initial_sha256"]
            or digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]):
        raise ValueError("Changed or unaccepted native state")
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
            raise ValueError("Wrong static active set")
        local_weight = saved["weight"].copy()
    with np.load(args.volume_rates) as saved:
        if not np.array_equal(saved["active"], active):
            raise ValueError("Wrong volume-rate active set")
        local_volume_rate = saved["weight_rate"].copy()
    weight = np.zeros(N)
    weight[active] = local_weight
    volume_rate = np.zeros(N)
    volume_rate[active] = local_volume_rate
    current_area, geometry, triangles, centers = power_areas(
        xy, weight, active, lookup, return_geometry=True)
    if geometry["hidden_active"] or np.any(current_area <= 0):
        raise ValueError("Invalid static power diagram")
    i, j, c, length = finite_faces(triangles, centers, lookup >= 0)
    d = np.linalg.norm(xy[j]-xy[i], axis=1)
    normal = (xy[j]-xy[i])/d[:, None]
    fixed_speed = (np.einsum("ij,ij->i", xy[j], velocity[j])
                   - np.einsum("ij,ij->i", xy[i], velocity[i])
                   - np.einsum("ij,ij->i", velocity[j]-velocity[i], c))/d
    v_kepler = np.column_stack((-c[:, 1], c[:, 0]))
    v_kepler *= (np.einsum("ij,ij->i", c, c)+EPS*EPS)[:, None]**(-.75)
    material_speed = np.einsum("ij,ij->i", v_kepler, normal)
    fixed_slip = fixed_speed-material_speed
    a = 1./(2*d)
    li, lj = lookup[i], lookup[j]
    conductance = length*a*a
    diag = np.zeros(len(active))
    rhs = np.zeros(len(active))
    np.add.at(diag, li, conductance)
    np.add.at(diag, lj, conductance)
    np.add.at(rhs, li, -length*a*fixed_slip)
    np.add.at(rhs, lj, length*a*fixed_slip)
    rows = np.r_[li, lj]
    cols = np.r_[lj, li]
    values = np.r_[-conductance, -conductance]
    lap = coo_matrix((values, (rows, cols)), shape=(len(active), len(active))).tocsr()
    lap.setdiag(diag)
    reduced = lap[1:, 1:]
    iterations = [0]
    solution, info = cg(
        reduced, rhs[1:],
        M=LinearOperator(reduced.shape, matvec=lambda v: v/diag[1:]),
        rtol=1e-9, atol=0., maxiter=3000,
        callback=lambda _: iterations.__setitem__(0, iterations[0]+1))
    if info != 0:
        raise RuntimeError(f"Material-rate CG failed after {iterations[0]} iterations: {info}")
    local_material_rate = np.r_[0., solution]
    local_material_rate -= np.mean(local_material_rate)
    material_rate = np.zeros(N)
    material_rate[active] = local_material_rate
    volume_slip = fixed_slip+a*(volume_rate[i]-volume_rate[j])
    material_slip = fixed_slip+a*(material_rate[i]-material_rate[j])
    same = np.abs(r0[i]-r0[j]) < 1e-10
    fixed_vdot, volume_vdot, _ = face_volume_rate(
        xy, velocity, weight, volume_rate, active, lookup, triangles, centers)
    _, material_vdot, _ = face_volume_rate(
        xy, velocity, weight, material_rate, active, lookup, triangles, centers)
    next_xy = xy+STEP*velocity
    target = p0["volume"][active].astype(float)
    predictions = {}
    for label, rate in (("volume_rate", volume_rate), ("material_rate", material_rate)):
        predicted_weight = weight+STEP*rate
        predicted_area, predicted_meta = power_areas(next_xy, predicted_weight, active, lookup)
        if predicted_meta["hidden_active"] or np.any(predicted_area <= 0):
            predictions[label] = {"valid_power_cells": False}
        else:
            predictions[label] = {"valid_power_cells": True,
                                  "area_change_from_current": area_metrics(predicted_area, current_area, r0[active]),
                                  "area_error_to_initial_target": area_metrics(predicted_area, target, r0[active])}
    result = {"accepted_receipt": str(args.accepted_receipt.resolve()),
              "endpoint_sha256": receipt["endpoint_sha256"],
              "static_weights_sha256": digest(args.static_weights),
              "volume_rates_sha256": digest(args.volume_rates),
              "time": 6.5, "disk_disk_faces": len(i),
              "objective": "Minimize sum over disk-disk faces of L*(v_face,n-v_Kepler(c).n)^2 with one scalar Wdot per disk cell",
              "material_rate_cg_iterations": iterations[0],
              "material_rate_linear_residual": float(np.linalg.norm(lap@local_material_rate-rhs)/np.linalg.norm(rhs)),
              "material_rate_abs_quantiles": np.quantile(np.abs(local_material_rate), [.5,.95,.99,1.]).tolist(),
              "fixed_weight_slip": grouped(length, fixed_slip, same),
              "volume_rate_slip": grouped(length, volume_slip, same),
              "material_rate_slip": grouped(length, material_slip, same),
              "volume_rate_cell_volume_rate_l2_over_fixed": float(np.linalg.norm(volume_vdot)/np.linalg.norm(fixed_vdot)),
              "material_rate_cell_volume_rate_l2_over_fixed": float(np.linalg.norm(material_vdot)/np.linalg.norm(fixed_vdot)),
              "euler_step": STEP, "one_step_geometry": predictions,
              "scope": "One geometry-only least-squares weight-rate comparison. Analytic Kepler material velocity, disk-disk faces only; no native weighted flux or long-time physical solution."}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    with args.rates.open("xb") as stream:
        np.savez_compressed(stream, active=active, material_rate=local_material_rate)
    for key in ("fixed_weight_slip", "volume_rate_slip", "material_rate_slip"):
        print(key, result[key]["all_disk_disk"]["length_weighted_mean_absolute_slip"])
    print("material cell Vdot / fixed", result["material_rate_cell_volume_rate_l2_over_fixed"],
          "CG", iterations[0])


if __name__ == "__main__":
    main()

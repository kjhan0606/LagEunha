"""Reconstruct weighted local volume rate and pressure-neutral correction."""

import argparse
import json
from pathlib import Path
import re

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import LinearOperator, cg

from audit_kepler_fresh_weighted_flux_proxy_t4952_20261027 import touching_finite_faces
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


BINS = ((.5, .52), (.52, .55), (.55, .65), (.5, .65),
        (.65, .75), (1., 2.), (.5, 2.))


def scatter_pair(i, j, flux, n):
    result = np.zeros(n)
    np.add.at(result, i, flux)
    np.add.at(result, j, -flux)
    return result


def compare(x, y, mass):
    w = mass/mass.sum()
    xx = float(np.sum(w*x*x))
    yy = float(np.sum(w*y*y))
    xy = float(np.sum(w*x*y))
    return {"predicted_rms": xx**.5, "observed_rms": yy**.5,
            "difference_rms": float(np.sqrt(np.sum(w*(x-y)**2))),
            "weighted_cosine": xy/(xx*yy)**.5 if xx and yy else None,
            "observed_on_predicted_slope": xy/xx if xx else None}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    audit = json.loads(args.receipt.read_text())
    if (len(audit["checks"]) != 21 or not all(audit["checks"].values())
            or audit["parent_time"] != 6.5
            or not 0. < audit["time"]-audit["parent_time"] <= .001000001):
        raise ValueError("Expected accepted native one-step receipt")
    for key, hash_key in (("initial", "initial_sha256"),
                          ("parent", "parent_sha256"),
                          ("endpoint", "endpoint_sha256")):
        if digest(Path(audit[key])) != audit[hash_key]:
            raise ValueError(f"Changed {key} snapshot")
    n = audit["n"]
    native_run = Path(audit["parent"]).parent
    p0, _ = ordered(Path(audit["initial"]), n)
    start, _ = ordered(Path(audit["parent"]), n)
    end, _ = ordered(Path(audit["endpoint"]), n)
    dt = float(audit["time"]-audit["parent_time"])
    active_mask = np.fromfile(native_run/"disk_active_mask.bin", dtype=np.uint8)
    if len(active_mask) != n or int(np.count_nonzero(active_mask)) != 49038:
        raise ValueError("Unexpected active disk mask")
    active = np.flatnonzero(active_mask)
    lookup = np.full(n, -1, dtype=int)
    lookup[active] = np.arange(len(active))
    r0 = np.hypot(p0["x"].astype(float)-3., p0["y"].astype(float)-3.)
    mass = p0["mass"].astype(float)
    xy = np.column_stack((start["x"].astype(float)-3.,
                          start["y"].astype(float)-3.))
    velocity = np.column_stack((start["vx"].astype(float),
                                start["vy"].astype(float)))
    pressure = start["pressure"].astype(float)
    weight = start["w2"].astype(float)
    area, meta, triangles, centers = power_areas(
        xy, weight, active, lookup, return_geometry=True)
    if meta["hidden_active"]:
        raise ValueError("Hidden active cell")
    area_relative_max = float(np.max(np.abs(
        area/start["volume"][active].astype(float)-1.)))
    if area_relative_max > 1e-6:
        raise ValueError("Independent power areas differ from native volumes")
    i, j, midpoint, length = touching_finite_faces(triangles, centers, active_mask != 0)
    dr = xy[j]-xy[i]
    d = np.linalg.norm(dr, axis=1)
    g = length/(2*d)
    h = scatter_pair(i, j, g*(pressure[i]-pressure[j]), n)
    h[active_mask == 0] = 0.
    jh = scatter_pair(i, j, g*(h[i]-h[j]), n)

    log = (native_run/"log").read_text(errors="replace")
    alpha = [float(x) for x in re.findall(
        r"\[GPU_VOLW2\] stage=[0-3] active=49038 faces=\d+ alpha=([0-9.eE+-]+)", log)]
    if len(alpha) != 4:
        raise ValueError("Expected four native GPU pressure-neutral stage projections")
    projected = -np.mean(alpha)*jh*dt/start["volume"].astype(float)
    observed = end["volume"].astype(float)/start["volume"].astype(float)-1.

    fixed_speed = (np.einsum("ij,ij->i", xy[j], velocity[j])
                   - np.einsum("ij,ij->i", xy[i], velocity[i])
                   - np.einsum("ij,ij->i", velocity[j]-velocity[i], midpoint))/d
    fixed_vdot = scatter_pair(i, j, length*fixed_speed, n)
    diagonal = np.zeros(n)
    np.add.at(diagonal, i, g)
    np.add.at(diagonal, j, g)
    both = (lookup[i] >= 0) & (lookup[j] >= 0)
    rows = np.r_[lookup[i[both]], lookup[j[both]]]
    columns = np.r_[lookup[j[both]], lookup[i[both]]]
    values = np.r_[-g[both], -g[both]]
    laplacian = coo_matrix((values, (rows, columns)),
                           shape=(len(active), len(active))).tocsr()
    laplacian += diags(diagonal[active], format="csr")
    rhs = -fixed_vdot[active]
    raw_rate = np.zeros(len(active))
    for _ in range(64):
        raw_rate += .8*(rhs-laplacian@raw_rate)/diagonal[active]
    h_active = h[active]
    alpha_cpu = float(np.dot(h_active, raw_rate)/np.dot(h_active, h_active))
    raw_residual = fixed_vdot[active]+laplacian@raw_rate
    projected_residual = raw_residual-alpha_cpu*(laplacian@h_active)
    correction_jacobi = np.zeros(len(active))
    for _ in range(64):
        correction_jacobi += .8*(h_active-laplacian@correction_jacobi)/diagonal[active]
    correction_exact, cg_info = cg(
        laplacian, h_active, rtol=1e-10, atol=0., maxiter=4000,
        M=LinearOperator(laplacian.shape,
                         matvec=lambda v: v/diagonal[active]))
    if cg_info != 0:
        raise ValueError(f"Inverse-gradient correction CG failed: {cg_info}")
    raw_work_rate = float(np.dot(h_active, raw_rate))
    denominator64 = float(np.dot(h_active, correction_jacobi))
    denominator_exact = float(np.dot(h_active, correction_exact))
    beta64 = raw_work_rate/denominator64
    beta_exact = raw_work_rate/denominator_exact
    inverse_jacobi_residual = raw_residual-beta64*(laplacian@correction_jacobi)
    inverse_exact_residual = raw_residual-beta_exact*h_active
    predicted_raw = np.zeros(n)
    predicted_cpu_projected = np.zeros(n)
    predicted_inverse_jacobi = np.zeros(n)
    predicted_inverse_exact = np.zeros(n)
    predicted_raw[active] = dt*raw_residual/start["volume"][active].astype(float)
    predicted_cpu_projected[active] = (
        dt*projected_residual/start["volume"][active].astype(float))
    predicted_inverse_jacobi[active] = (
        dt*inverse_jacobi_residual/start["volume"][active].astype(float))
    predicted_inverse_exact[active] = (
        dt*inverse_exact_residual/start["volume"][active].astype(float))
    secant_rate = (end["w2"].astype(float)-weight)/dt
    secant_rate[active_mask == 0] = 0.
    secant_vdot = scatter_pair(i, j, g*(secant_rate[i]-secant_rate[j]), n)
    kinematic = dt*(fixed_vdot+secant_vdot)/start["volume"].astype(float)

    bins = {}
    for lo, hi in BINS:
        take = (r0 >= lo) & (r0 < hi)
        w = mass[take]
        bins[f"{lo:g}:{hi:g}"] = {
            "cells": int(take.sum()),
            "observed_mean_abs_volume_change": float(np.sum(w*np.abs(observed[take]))/w.sum()),
            "observed_p99_abs_volume_change": float(np.quantile(np.abs(observed[take]), .99)),
            "pressure_neutral_projection_vs_observed": compare(projected[take], observed[take], w),
            "cpu_jacobi64_raw_vs_observed": compare(predicted_raw[take], observed[take], w),
            "cpu_jacobi64_projected_vs_observed": compare(
                predicted_cpu_projected[take], observed[take], w),
            "inverse_gradient_jacobi64_predicted_rms": compare(
                predicted_inverse_jacobi[take], observed[take], w)["predicted_rms"],
            "inverse_gradient_exact_predicted_rms": compare(
                predicted_inverse_exact[take], observed[take], w)["predicted_rms"],
            "secant_kinematic_vs_observed": compare(kinematic[take], observed[take], w),
        }
    report = {"receipt": str(args.receipt), "receipt_sha256": digest(args.receipt),
              "start_sha256": audit["parent_sha256"],
              "endpoint_sha256": audit["endpoint_sha256"],
              "time": audit["parent_time"], "dt": dt,
              "independent_power_area_max_relative_error": area_relative_max,
              "finite_faces_touching_disk": len(i),
              "stage_projection_alpha": alpha,
              "cpu_jacobi64_alpha": alpha_cpu,
              "cpu_alpha_relative_gap_to_native_stage0": alpha_cpu/alpha[0]-1.,
              "inverse_gradient_jacobi64_beta": beta64,
              "inverse_gradient_exact_beta": beta_exact,
              "inverse_gradient_jacobi64_relative_linear_residual": float(
                  np.linalg.norm(laplacian@correction_jacobi-h_active)
                  /np.linalg.norm(h_active)),
              "inverse_gradient_exact_relative_linear_residual": float(
                  np.linalg.norm(laplacian@correction_exact-h_active)
                  /np.linalg.norm(h_active)),
              "inverse_gradient_work_neutrality_jacobi64": float(
                  np.dot(h_active, raw_rate-beta64*correction_jacobi)),
              "inverse_gradient_work_neutrality_exact": float(
                  np.dot(h_active, raw_rate-beta_exact*correction_exact)),
              "unprojected_pressure_weight_work_one_step": float(
                  dt*np.dot(h_active, raw_rate)),
              "unprojected_pressure_weight_work_over_internal_energy": float(
                  dt*np.dot(h_active, raw_rate)/np.sum(start["ie"].astype(float))),
              "bins": bins,
              "scope": ("Independent t=6.5 regular triangulation and saved "
                        "one-step weight secant. Projection-only prediction "
                        "uses the mean logged RK4 alpha times J h on the "
                        "initial geometry. CPU Jacobi64 reproduces the source "
                        "iteration and scalar projection on independent "
                        "geometry. Both omit RK4 stage variation and are "
                        "not a captured native stage volume rate. Secant "
                        "kinematics are a separate first-order check with "
                        "potential strong cancellation errors. Inverse-gradient "
                        "projection is an unevolved constrained-rate candidate, "
                        "not an observed native evolution result.")}
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    for key in ("0.5:0.52", "0.52:0.55", "0.55:0.65", "0.5:2"):
        row = bins[key]
        print(key, "obs", row["observed_mean_abs_volume_change"],
              "proj", row["pressure_neutral_projection_vs_observed"],
              "jacobi_raw", row["cpu_jacobi64_raw_vs_observed"],
              "jacobi_projected", row["cpu_jacobi64_projected_vs_observed"],
              "inverse64_rms", row["inverse_gradient_jacobi64_predicted_rms"],
              "inverse_exact_rms", row["inverse_gradient_exact_predicted_rms"],
              "secant", row["secant_kinematic_vs_observed"])
    print("native_alpha0", alpha[0], "cpu_alpha", alpha_cpu,
          "unprojected_work_over_internal_energy",
          report["unprojected_pressure_weight_work_over_internal_energy"])


if __name__ == "__main__":
    main()

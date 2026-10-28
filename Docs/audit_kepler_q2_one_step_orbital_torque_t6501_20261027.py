"""Measure native one-step orbital angular-momentum change by initial ring."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


BINS = ((.5, .52), (.52, .55), (.55, .65), (.5, .65),
        (.65, .75), (.75, 1.), (1., 2.), (.5, 2.))
EPS = .01175


def specific_orbital_l(p):
    return ((p["x"].astype(float)-3.)*p["vy"].astype(float)
            -(p["y"].astype(float)-3.)*p["vx"].astype(float))


def weighted_mean(v, w):
    return float(np.sum(v*w, dtype=np.longdouble)/np.sum(w, dtype=np.longdouble))


def gravity_only_rk4_l_change(p, dt, ell0):
    """Independent pure-central-gravity RK4 reference from a saved state."""
    r = np.column_stack((p["x"].astype(float)-3., p["y"].astype(float)-3.))
    v = np.column_stack((p["vx"].astype(float), p["vy"].astype(float)))

    def gravity(pos):
        radius2 = np.einsum("ij,ij->i", pos, pos)
        return -pos*(radius2+EPS*EPS)[:, None]**(-1.5)

    k1r, k1v = v, gravity(r)
    k2r = v+.5*dt*k1v
    k2v = gravity(r+.5*dt*k1r)
    k3r = v+.5*dt*k2v
    k3v = gravity(r+.5*dt*k2r)
    k4r = v+dt*k3v
    k4v = gravity(r+dt*k3r)
    end_r = r+dt*(k1r+2*k2r+2*k3r+k4r)/6.
    end_v = v+dt*(k1v+2*k2v+2*k3v+k4v)/6.
    end_l = end_r[:, 0]*end_v[:, 1]-end_r[:, 1]*end_v[:, 0]
    return (end_l-specific_orbital_l(p))/ell0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("w0_t65_receipt", type=Path)
    parser.add_argument("w0_step_receipt", type=Path)
    parser.add_argument("weighted_t65_receipt", type=Path)
    parser.add_argument("weighted_step_receipt", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    w0_t65, w0_step, weighted_t65, weighted_step = [json.loads(path.read_text())
        for path in (args.w0_t65_receipt, args.w0_step_receipt,
                     args.weighted_t65_receipt, args.weighted_step_receipt)]
    if not all(w0_t65["checks"].values()):
        raise ValueError("W0 t6.5 gate failed")
    if not all(w0_step["checks"].values()):
        raise ValueError("W0 one-step gate failed")
    if not all(weighted_t65["checks"].values()):
        raise ValueError("Weighted t6.5 gate failed")
    if not all(weighted_step["checks"].values()):
        raise ValueError("Weighted one-step gate failed")
    if (w0_t65["initial_sha256"] != weighted_t65["initial_sha256"]
            or w0_t65["initial_sha256"] != weighted_step["initial_sha256"]):
        raise ValueError("Different initial Kepler IC")
    if (w0_t65["endpoint_sha256"] != w0_step["restart_sha256"]
            or weighted_t65["endpoint_sha256"] != weighted_step["parent_sha256"]):
        raise ValueError("One-step start does not match accepted t6.5 endpoint")
    if (w0_step["start_time"] != 6.5 or w0_step["endpoint_time"] != 6.501
            or weighted_step["parent_time"] != 6.5
            or weighted_step["time"] != 6.501):
        raise ValueError("Expected identical t6.5 to t6.501 intervals")

    initial_path = Path(w0_t65["initial"])
    if digest(initial_path) != w0_t65["initial_sha256"]:
        raise ValueError("Initial snapshot digest mismatch")
    n = w0_t65["n"]
    p0, _ = ordered(initial_path, n)
    r0 = np.hypot(p0["x"].astype(float)-3., p0["y"].astype(float)-3.)
    mass = p0["mass"].astype(float)
    ell0 = specific_orbital_l(p0)
    result = {"initial_sha256": w0_t65["initial_sha256"],
              "start_time": 6.5, "endpoint_time": 6.501,
              "dt": .001, "bins": {}}
    cases = {
        "w0": (w0_step["restart"], w0_step["restart_sha256"],
               w0_step["endpoint"], w0_step["endpoint_sha256"]),
        "weighted": (weighted_step["parent"], weighted_step["parent_sha256"],
                     weighted_step["endpoint"], weighted_step["endpoint_sha256"]),
    }
    for label, (start_path, start_sha, end_path, end_sha) in cases.items():
        start_path, end_path = Path(start_path), Path(end_path)
        if digest(start_path) != start_sha or digest(end_path) != end_sha:
            raise ValueError(f"{label}: changed native snapshot")
        a, ma = ordered(start_path, n)
        b, mb = ordered(end_path, n)
        if (abs(ma["t"]-6.5) > 1e-12 or abs(mb["t"]-6.501) > 1e-12
                or not np.array_equal(a["mass"], p0["mass"])
                or not np.array_equal(b["mass"], p0["mass"])):
            raise ValueError(f"{label}: time or mass mismatch")
        start_error = specific_orbital_l(a)/ell0-1.
        end_error = specific_orbital_l(b)/ell0-1.
        increment = end_error-start_error
        gravity_reference = gravity_only_rk4_l_change(a, .001, ell0)
        bins = {}
        for lo, hi in BINS:
            take = (r0 >= lo) & (r0 < hi)
            w = mass[take]
            de = increment[take]
            s = start_error[take]
            e = end_error[take]
            bins[f"{lo:g}:{hi:g}"] = {
                "cells": int(take.sum()),
                "mean_abs_increment_relative_to_initial_l": weighted_mean(np.abs(de), w),
                "p99_abs_increment_relative_to_initial_l": float(np.quantile(np.abs(de), .99)),
                "signed_mean_increment_relative_to_initial_l": weighted_mean(de, w),
                "mean_abs_error_at_start": weighted_mean(np.abs(s), w),
                "mean_abs_error_at_end": weighted_mean(np.abs(e), w),
                "mean_change_in_abs_error": weighted_mean(np.abs(e)-np.abs(s), w),
                "same_sign_as_prior_error_mass_fraction": weighted_mean(
                    ((de*s) > 0).astype(float), w),
                "pure_gravity_RK4_mean_abs_increment": weighted_mean(
                    np.abs(gravity_reference[take]), w),
                "pure_gravity_RK4_p99_abs_increment": float(np.quantile(
                    np.abs(gravity_reference[take]), .99)),
            }
        result["bins"][label] = bins
        result[label+"_start_sha256"] = start_sha
        result[label+"_endpoint_sha256"] = end_sha
    result["scope"] = ("Native one-step orbital angular-momentum increment "
                       "from accepted snapshots. Central analytic gravity "
                       "has zero continuous torque, but the increment also "
                       "contains RK4 and snapshot rounding effects. An "
                       "independent pure-central-gravity RK4 one-step "
                       "reference is reported, not a native force-off run. W0 and "
                       "weighted cases have the same t0 IC but different "
                       "evolved t6.5 states; this is not a same-state force ablation.")
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    for key in ("0.5:0.52", "0.52:0.55", "0.55:0.65", "0.5:2"):
        print(key, "W0", result["bins"]["w0"][key]["mean_abs_increment_relative_to_initial_l"],
              "weighted", result["bins"]["weighted"][key]["mean_abs_increment_relative_to_initial_l"])


if __name__ == "__main__":
    main()

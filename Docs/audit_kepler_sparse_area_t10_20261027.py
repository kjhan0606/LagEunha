"""Audit corrected sparse q2 native Kepler at an exact time (default t=10)."""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered, mean_abs
from prepare_kepler_sparse_area_restart_20261027 import IC, digest, snapshot_time


T = 10.
EPS = .01175
GAMMA_MINUS_ONE = float(np.float32(1.666667)) - 1.
BINS = ((0., .5), (.5, .75), (.75, 1.), (1., 2.), (.5, 2.),
        (2., 2.35), (2.35, 4.25))


def energy(p):
    x = p["x"].astype(float) - 3.
    y = p["y"].astype(float) - 3.
    m = p["mass"].astype(float)
    vx, vy = p["vx"].astype(float), p["vy"].astype(float)
    kinetic = np.sum(.5 * m * (vx * vx + vy * vy), dtype=np.longdouble)
    thermal = np.sum(p["ie"].astype(float), dtype=np.longdouble)
    potential = np.sum(-m / np.sqrt(x*x + y*y + EPS*EPS), dtype=np.longdouble)
    return float(kinetic + thermal + potential)


def ledger_last(run):
    # Long runs can write hundreds of MB of tessellation diagnostics. The
    # final RK4 ledger is near EOF; avoid retaining the entire log in RAM.
    with (run / "log").open("rb") as stream:
        stream.seek(0, 2)
        stream.seek(max(0, stream.tell() - 8 * 1024 * 1024))
        tail = stream.read().decode("utf-8", errors="replace")
    lines = [line for line in tail.splitlines() if "[RK4E]" in line]
    if not lines:
        raise ValueError("Missing native RK4 energy ledger")
    result = {}
    for key in ("floor_cum", "de_cum", "sfl_cum"):
        found = re.search(rf"\b{key}=([0-9.eE+-]+)", lines[-1])
        if not found:
            raise ValueError(f"Missing final RK4 {key}")
        result[key] = float(found.group(1))
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("initial", type=Path)
    ap.add_argument("run", type=Path)
    ap.add_argument("--expected-count", type=int, choices=tuple(IC), required=True)
    ap.add_argument("--time", type=float, default=T)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    target = args.time
    if not np.isfinite(target) or target <= 0:
        raise ValueError("Expected finite positive target time")
    n = args.expected_count
    if abs(snapshot_time(args.initial, n)) > 1e-12:
        raise ValueError("Reference is not an exact t=0 native checkpoint")
    candidates = [path for path in args.run.glob("khout.*.dat")
                  if abs(snapshot_time(path, n) - target) < 1e-10]
    if len(candidates) != 1:
        raise ValueError("Missing unique complete exact target-time checkpoint")
    endpoint = candidates[0]
    p0, m0 = ordered(args.initial, n)
    p, mt = ordered(endpoint, n)
    if abs(m0["t"]) > 1e-12 or abs(mt["t"] - target) > 1e-10:
        raise ValueError("Native snapshot time mismatch")

    x0, y0 = p0["x"].astype(float)-3., p0["y"].astype(float)-3.
    x, y = p["x"].astype(float)-3., p["y"].astype(float)-3.
    r0, r = np.hypot(x0, y0), np.hypot(x, y)
    m = p0["mass"].astype(float)
    phase = target * (r0*r0 + EPS*EPS)**(-.75)
    xb = x0*np.cos(phase) - y0*np.sin(phase)
    yb = x0*np.sin(phase) + y0*np.cos(phase)
    position_error = np.hypot((x-xb+3.)%6.-3., (y-yb+3.)%6.-3.)
    density_change = p["den"].astype(float) / p0["den"].astype(float) - 1.
    pressure_change = p["pressure"].astype(float) / p0["pressure"].astype(float) - 1.
    entropy_log = (np.log(p["ie"].astype(float) / p0["ie"].astype(float)) +
                   GAMMA_MINUS_ONE * np.log(p["volume"].astype(float) /
                                            p0["volume"].astype(float)))
    vr = (x*p["vx"].astype(float) + y*p["vy"].astype(float)) / np.maximum(r, 1e-300)
    vk0 = r0 * (r0*r0 + EPS*EPS)**(-.75)
    radial_bins = []
    for lo, hi in BINS:
        take = (r0 >= lo) & (r0 < hi)
        good = take & (r0 >= .05)
        if not np.any(take):
            continue
        weights = m[take]
        l0 = np.sum(weights*(x0[take]*p0["vy"][take] - y0[take]*p0["vx"][take]),
                    dtype=np.longdouble)
        lt = np.sum(weights*(x[take]*p["vy"][take] - y[take]*p["vx"][take]),
                    dtype=np.longdouble)
        radial_bins.append({
            "rlo": lo, "rhi": hi, "cells": int(np.sum(take)),
            "mass": float(np.sum(weights, dtype=np.longdouble)),
            "density_abs_mass_weighted": mean_abs(density_change[take], weights),
            "pressure_abs_mass_weighted": mean_abs(pressure_change[take], weights),
            "entropy_log_signed_mass_weighted": float(np.sum(entropy_log[take]*weights,
                                                            dtype=np.longdouble)/np.sum(weights, dtype=np.longdouble)),
            "entropy_log_abs_mass_weighted": mean_abs(entropy_log[take], weights),
            "position_error_abs_mass_weighted": mean_abs(position_error[take], weights),
            "position_error_p95": float(np.quantile(position_error[take], .95)),
            "radial_speed_rms_over_circular": float(np.sqrt(
                np.sum((vr[good]/vk0[good])**2*m[good], dtype=np.longdouble) /
                np.sum(m[good], dtype=np.longdouble))) if np.any(good) else None,
            "angular_momentum_relative": float(lt/l0-1.) if l0 else None,
        })
    disk = next(row for row in radial_bins if row["rlo"] == .5 and row["rhi"] == 2.)
    disk_mask = (r0 >= .5) & (r0 < 2.)
    annulus = (r > .8) & (r < 1.2)
    if not np.any(annulus):
        raise ValueError("Empty paper-reference annulus")
    vphi = (x[annulus]*p["vy"][annulus] - y[annulus]*p["vx"][annulus])/r[annulus]
    vka = r[annulus]*(r[annulus]**2 + EPS**2)**(-.75)
    l1v = float(np.mean(np.abs(vphi-vka)/vka))
    l1rho_literal = float(np.mean(np.abs(p["den"][annulus].astype(float)-1.)))
    l1rho_plateau = float(np.mean(np.abs(p["den"][annulus].astype(float)-1.01)/1.01))
    e0, et = energy(p0), energy(p)
    edrift = (et-e0)/abs(e0)
    ledger = ledger_last(args.run)
    volume0 = float(np.sum(p0["volume"].astype(float), dtype=np.longdouble))
    volumet = float(np.sum(p["volume"].astype(float), dtype=np.longdouble))
    gm1 = GAMMA_MINUS_ONE
    rho_max_relative = max(float(np.max(np.abs(
        q["den"].astype(float) /
        (q["mass"].astype(float) / q["volume"].astype(float)) - 1.)))
        for q in (p0, p))
    eos_pressure_max_relative = max(float(np.max(np.abs(
        q["pressure"].astype(float) /
        (gm1 * q["ie"].astype(float) / q["volume"].astype(float)) - 1.)))
        for q in (p0, p))
    finite = all(np.isfinite(q[field]).all() for q in (p0, p) for field in
                 ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2"))
    checks = {
        "finite": bool(finite),
        "initial_periodic_volume": abs(volume0-36.) < 1e-6,
        "endpoint_periodic_volume": abs(volumet-36.) < 1e-5,
        "mass_exact": bool(np.array_equal(p0["mass"], p["mass"])),
        "W0_exact": bool(np.all(p0["w2"] == 0) and np.all(p["w2"] == 0)),
        "rho_equals_mass_over_volume": bool(np.isfinite(rho_max_relative) and
                                             rho_max_relative < 1e-6),
        "pressure_matches_ideal_gas_EOS": bool(np.isfinite(eos_pressure_max_relative) and
                                                 eos_pressure_max_relative < 1e-6),
        "positive_state": bool(np.all(p0["volume"] > 0) and np.all(p["volume"] > 0) and
                               np.all(p0["ie"] > 0) and np.all(p["ie"] > 0)),
        "no_floor_injection": all(value == 0 for value in ledger.values()),
        "total_energy_relative_below_1e_minus_3": bool(np.isfinite(edrift) and abs(edrift) < .001),
        "gizmo_annulus_L1_velocity_below_0p01": l1v < .01,
        "gizmo_annulus_L1_density_below_0p25": l1rho_literal < .25,
        "disk_density_abs_below_0p05": disk["density_abs_mass_weighted"] < .05,
        "disk_entropy_log_abs_below_0p01": disk["entropy_log_abs_mass_weighted"] < .01,
        "disk_angular_momentum_drift_below_0p01": abs(disk["angular_momentum_relative"]) < .01,
    }
    result = {
        "initial": str(args.initial), "endpoint": str(endpoint),
        "endpoint_sha256": digest(endpoint), "time": float(mt["t"]), "n": n,
        "checks": checks, "initial_volume_sum": volume0, "endpoint_volume_sum": volumet,
        "mass_sum": float(np.sum(m, dtype=np.longdouble)),
        "disk_cells": int(np.sum(disk_mask)), "gizmo_annulus_cells": int(np.sum(annulus)),
        "gizmo_L1_velocity": l1v,
        "gizmo_L1_density_literal": l1rho_literal,
        "gizmo_L1_density_plateau_corrected": l1rho_plateau,
        "total_energy_relative": edrift, "initial_energy": e0, "endpoint_energy": et,
        "rho_mass_volume_max_relative": rho_max_relative,
        "eos_pressure_max_relative": eos_pressure_max_relative,
        "rk4_ledger": ledger, "minimum_internal_energy": float(np.min(p["ie"])),
        "minimum_volume": float(np.min(p["volume"])), "radial_bins": radial_bins,
        "scope": f"Corrected sparse W0 native CUDA+MPI+OpenMP t={target:g} gate; not proof of later stability",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))
    if not all(checks.values()):
        raise SystemExit("Corrected sparse t=10 physical gate failed")


if __name__ == "__main__":
    main()

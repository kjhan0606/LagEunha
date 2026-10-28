"""Gate a complete native polar Kepler restart at an exact later time."""

import argparse
import json
from pathlib import Path
import re
import struct
import sys

import numpy as np

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot


N = 384 * 384
EPS = .01175
CX = CY = 3.
GM1 = float(np.float32(1.666667)) - 1.
BINS = ((0., .1), (.1, .5), (.5, .75), (.75, 1.), (1., 2.),
        (.5, 2.), (2., 2.35), (2.35, 2.75), (2.75, 5.))


def snapshot_time(path):
    with path.open("rb") as stream:
        magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
        stream.seek(-52, 2)
        t = struct.unpack("<d", stream.read(8))[0]
    if ((magic, n, stride, stress) != (-1, N, 576, 1)
            or path.stat().st_size != 16 + n * stride + 52):
        raise ValueError(f"Incomplete or incompatible native checkpoint: {path}")
    return t


def ordered(path):
    time = snapshot_time(path)
    p, meta = read_snapshot(path)
    ids = (p["u4if"].astype(np.uint64) & np.uint64((1 << 56) - 1)).astype(np.int64)
    if len(p) != N or not np.array_equal(np.sort(ids), np.arange(N)):
        raise ValueError(f"Incomplete cell ID set: {path}")
    if abs(float(meta["t"]) - time) > 1e-12:
        raise ValueError("Header/footer time disagreement")
    return p[np.argsort(ids)], time


def weighted_mean(value, weight):
    return float(np.sum(value * weight, dtype=np.longdouble)
                 / np.sum(weight, dtype=np.longdouble))


def total_energy(p):
    x = p["x"].astype(float) - CX
    y = p["y"].astype(float) - CY
    m = p["mass"].astype(float)
    vx = p["vx"].astype(float)
    vy = p["vy"].astype(float)
    kinetic = .5 * np.sum(m * (vx * vx + vy * vy), dtype=np.longdouble)
    thermal = np.sum(p["ie"].astype(float), dtype=np.longdouble)
    potential = -np.sum(m / np.sqrt(x * x + y * y + EPS * EPS),
                        dtype=np.longdouble)
    return float(kinetic + thermal + potential)


def final_ledger(run):
    lines = [line for line in (run / "log").read_text(errors="replace").splitlines()
             if "[RK4E]" in line]
    if not lines:
        raise ValueError("Missing native RK4 energy ledger")
    last = lines[-1]
    result = {}
    for key in ("floor_cum", "de_cum", "sfl_cum"):
        match = re.search(rf"\b{key}=([0-9.eE+-]+)", last)
        if not match:
            raise ValueError(f"Missing {key} in final native RK4 ledger")
        result[key] = float(match.group(1))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("run", type=Path)
    parser.add_argument("--time", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    p0, t0 = ordered(args.initial)
    if t0 != 0 or args.time <= 0:
        raise ValueError("Expected exact t=0 reference and positive target")
    matches = []
    for path in args.run.glob("khout.*.dat"):
        if abs(snapshot_time(path) - args.time) < 1e-10:
            matches.append(path)
    if len(matches) != 1:
        raise ValueError("Missing unique complete exact-time checkpoint")
    endpoint = matches[0]
    p, t = ordered(endpoint)
    m = p0["mass"].astype(float)
    x0 = p0["x"].astype(float) - CX
    y0 = p0["y"].astype(float) - CY
    r0 = np.hypot(x0, y0)
    x = p["x"].astype(float) - CX
    y = p["y"].astype(float) - CY
    r = np.hypot(x, y)
    omega = (r0 * r0 + EPS * EPS)**(-.75)
    phase = t * omega
    xb = CX + x0 * np.cos(phase) - y0 * np.sin(phase)
    yb = CY + x0 * np.sin(phase) + y0 * np.cos(phase)
    position_error = np.hypot(p["x"] - xb, p["y"] - yb)
    density_change = p["den"].astype(float) / p0["den"].astype(float) - 1.
    entropy_log = (np.log(p["ie"].astype(float) / p0["ie"].astype(float))
                   + GM1 * np.log(p["volume"].astype(float)
                                  / p0["volume"].astype(float)))
    vr = (x * p["vx"] + y * p["vy"]) / np.maximum(r, 1e-300)
    disk = (r0 >= .5) & (r0 < 2.)
    annulus = (r > .8) & (r < 1.2)
    if not np.any(annulus):
        raise ValueError("Empty GIZMO annulus")
    vphi = (x[annulus] * p["vy"][annulus]
            - y[annulus] * p["vx"][annulus]) / r[annulus]
    vk = r[annulus] * (r[annulus]**2 + EPS**2)**(-.75)
    l1_v = float(np.mean(abs(vphi - vk) / vk))
    l1_rho = float(np.mean(abs(p["den"][annulus].astype(float) - 1.)))
    l1_rho_plateau = float(np.mean(abs(p["den"][annulus].astype(float) - 1.01) / 1.01))
    L0 = np.sum(m[disk] * (x0[disk] * p0["vy"][disk]
                           - y0[disk] * p0["vx"][disk]), dtype=np.longdouble)
    Lt = np.sum(m[disk] * (x[disk] * p["vy"][disk]
                           - y[disk] * p["vx"][disk]), dtype=np.longdouble)
    E0, Et = total_energy(p0), total_energy(p)
    ledger = final_ledger(args.run)
    area0 = float(np.sum(p0["volume"].astype(float), dtype=np.longdouble))
    areat = float(np.sum(p["volume"].astype(float), dtype=np.longdouble))
    rho_mass_volume_max_relative = max(
        float(np.max(abs(q["den"].astype(float) /
                         (q["mass"].astype(float) / q["volume"].astype(float)) - 1.)))
        for q in (p0, p))
    rows = []
    for lo, hi in BINS:
        mask = (r0 >= lo) & (r0 < hi)
        rows.append({"rlo": lo, "rhi": hi, "cells": int(np.sum(mask)),
                     "density_abs_mass_weighted": weighted_mean(abs(density_change[mask]), m[mask]),
                     "entropy_log_abs_mass_weighted": weighted_mean(abs(entropy_log[mask]), m[mask]),
                     "radial_speed_abs_mass_weighted": weighted_mean(abs(vr[mask]), m[mask]),
                     "position_error_abs_mass_weighted": weighted_mean(position_error[mask], m[mask])})
    disk_density = weighted_mean(abs(density_change[disk]), m[disk])
    disk_entropy = weighted_mean(abs(entropy_log[disk]), m[disk])
    disk_position = weighted_mean(position_error[disk], m[disk])
    Ldrift = float(Lt / L0 - 1.)
    Edrift = (Et - E0) / abs(E0)
    checks = {
        "finite": bool(all(np.isfinite(q[field]).all() for q in (p0, p) for field in
                           ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2"))),
        "mass_exact": bool(np.array_equal(p0["mass"], p["mass"])),
        "W0_exact": bool(np.all(p0["w2"] == 0) and np.all(p["w2"] == 0)),
        "positive_U_and_volume": bool(np.all(p0["ie"] > 0) and np.all(p["ie"] > 0)
                                      and np.all(p0["volume"] > 0) and np.all(p["volume"] > 0)),
        "initial_area_tiles_6x6_box": bool(abs(area0 - 36.) < 1e-5),
        "endpoint_area_tiles_6x6_box": bool(abs(areat - 36.) < 1e-5),
        "rho_equals_mass_over_volume": bool(np.isfinite(rho_mass_volume_max_relative)
                                             and rho_mass_volume_max_relative < 1e-6),
        "no_floor_or_dual_energy_injection": bool(all(value == 0 for value in ledger.values())),
        "gizmo_L1_velocity_below_0p01": bool(l1_v < .01),
        "gizmo_L1_density_below_0p25": bool(l1_rho < .25),
        "disk_density_change_below_0p05": bool(disk_density < .05),
        "disk_entropy_log_abs_below_0p01": bool(disk_entropy < .01),
        "disk_position_error_below_0p01": bool(disk_position < .01),
        "disk_angular_momentum_drift_below_0p01": bool(abs(Ldrift) < .01),
        "total_energy_drift_below_0p001": bool(abs(Edrift) < .001),
    }
    result = {
        "initial": str(args.initial), "endpoint": str(endpoint),
        "time": t, "n": N, "disk_cells": int(np.sum(disk)),
        "disk_density_abs_mass_weighted": disk_density,
        "disk_entropy_log_abs_mass_weighted": disk_entropy,
        "disk_entropy_log_signed_mass_weighted": weighted_mean(entropy_log[disk], m[disk]),
        "disk_position_error_abs_mass_weighted": disk_position,
        "disk_angular_momentum_relative": Ldrift,
        "total_energy_relative": Edrift,
        "initial_area_sum": area0,
        "endpoint_area_sum": areat,
        "rho_mass_over_volume_max_relative": rho_mass_volume_max_relative,
        "gizmo_annulus_cells": int(np.sum(annulus)),
        "gizmo_L1_velocity": l1_v,
        "gizmo_L1_density_literal": l1_rho,
        "gizmo_L1_density_plateau_corrected": l1_rho_plateau,
        "rk4_final_ledger": ledger,
        "minimum_internal_energy": float(np.min(p["ie"])),
        "minimum_volume": float(np.min(p["volume"])),
        "bins": rows, "checks": checks,
        "scope": "Native CUDA+MPI+OpenMP fixed-mass W=0 polar Kepler checkpoint; gate is not a t=600 stability claim",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))
    if not all(checks.values()):
        raise SystemExit("Native long-time physics gate failed")


if __name__ == "__main__":
    main()

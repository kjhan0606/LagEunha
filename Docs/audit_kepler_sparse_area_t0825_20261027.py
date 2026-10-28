"""Gate corrected sparse q2 native evolution at exact t=.825."""

import argparse
import json
from pathlib import Path
import re
import sys

import numpy as np

from prepare_kepler_sparse_area_restart_20261027 import IC, digest, snapshot_time

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot


TARGET = .825
GAMMA_MINUS_ONE = float(np.float32(1.666667)) - 1.
BINS = ((0., .2), (.2, .5), (0., .5), (.5, .75), (.75, 1.),
        (1., 2.), (.5, 2.), (2., 2.35), (2.35, 4.25))


def ordered(path, n):
    p, meta = read_snapshot(path)
    ids = (p["u4if"].astype(np.uint64) & np.uint64((1 << 56) - 1)).astype(np.int64)
    if len(p) != n or not np.array_equal(np.sort(ids), np.arange(n)):
        raise ValueError(f"Incomplete native IDs: {path}")
    return p[np.argsort(ids)], meta


def mean_abs(value, weights):
    return float(np.sum(np.abs(value) * weights, dtype=np.longdouble) /
                 np.sum(weights, dtype=np.longdouble))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("initial", type=Path)
    ap.add_argument("run", type=Path)
    ap.add_argument("--expected-count", type=int, choices=tuple(IC), required=True)
    ap.add_argument("--target-time", type=float, default=TARGET)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if not np.isfinite(args.target_time) or args.target_time <= 0:
        raise ValueError("Target time must be finite and positive")
    n = args.expected_count
    if abs(snapshot_time(args.initial, n)) > 1e-12:
        raise ValueError("Reference is not exact initial state")
    candidates = [path for path in args.run.glob("khout.*.dat")
                  if abs(snapshot_time(path, n) - args.target_time) < 1e-12]
    if len(candidates) != 1:
        raise ValueError("Missing unique complete exact target-time checkpoint")
    endpoint = candidates[0]
    p0, m0 = ordered(args.initial, n)
    p, m = ordered(endpoint, n)
    if abs(m0["t"]) > 1e-12 or abs(m["t"] - args.target_time) > 1e-12:
        raise ValueError("Native snapshot time mismatch")
    mass = p0["mass"].astype(float)
    r0x, r0y = p0["x"] - 3., p0["y"] - 3.
    r0 = np.hypot(r0x, r0y)
    phase = args.target_time * (r0 * r0 + .01175 ** 2) ** (-.75)
    orbit_x = 3 + r0x * np.cos(phase) - r0y * np.sin(phase)
    orbit_y = 3 + r0x * np.sin(phase) + r0y * np.cos(phase)
    position_error = np.hypot(p["x"] - orbit_x, p["y"] - orbit_y)
    rho_change = p["den"].astype(float) / p0["den"].astype(float) - 1
    pressure_change = p["pressure"].astype(float) / p0["pressure"].astype(float) - 1
    entropy_change = p["ie"].astype(float) / (
        p0["ie"].astype(float) *
        (p0["volume"].astype(float) / p["volume"].astype(float)) ** GAMMA_MINUS_ONE) - 1
    x1, y1 = p["x"] - 3., p["y"] - 3.
    radius = np.hypot(x1, y1)
    vr = (x1 * p["vx"] + y1 * p["vy"]) / np.maximum(radius, 1e-300)
    omega = (r0 * r0 + .01175 ** 2) ** (-.75)
    vcirc = r0 * omega
    bins = []
    for lo, hi in BINS:
        take = (r0 >= lo) & (r0 < hi)
        valid_v = take & (r0 >= .05)
        w = mass[take]
        wv = mass[valid_v]
        l0 = np.sum(w * (r0x[take] * p0["vy"][take] - r0y[take] * p0["vx"][take]), dtype=np.longdouble)
        l1 = np.sum(w * (x1[take] * p["vy"][take] - y1[take] * p["vx"][take]), dtype=np.longdouble)
        bins.append({"rlo": lo, "rhi": hi, "cells": int(np.sum(take)),
                     "mass": float(np.sum(w, dtype=np.longdouble)),
                     "density_abs_mass_weighted": mean_abs(rho_change[take], w),
                     "density_abs_p95": float(np.quantile(abs(rho_change[take]), .95)),
                     "pressure_abs_mass_weighted": mean_abs(pressure_change[take], w),
                     "entropy_signed_mass_weighted": float(np.sum(entropy_change[take] * w, dtype=np.longdouble) /
                                                           np.sum(w, dtype=np.longdouble)),
                     "entropy_abs_mass_weighted": mean_abs(entropy_change[take], w),
                     "position_error_mass_weighted": mean_abs(position_error[take], w),
                     "position_error_max": float(np.max(position_error[take])),
                     "radial_speed_rms_over_circular": float(np.sqrt(
                         np.sum((vr[valid_v] / vcirc[valid_v]) ** 2 * wv, dtype=np.longdouble) /
                         np.sum(wv, dtype=np.longdouble))),
                     "angular_momentum_relative": float(l1 / l0 - 1) if l0 != 0 else None})
    disk = next(row for row in bins if row["rlo"] == .5 and row["rhi"] == 2.)
    first_ring = next(row for row in bins if row["rlo"] == 0. and row["rhi"] == .2)
    inner0, innert = r0 < .5, radius < .5
    inward = ~inner0 & innert
    outward = inner0 & ~innert
    crossing_mass = float(np.sum(mass[inward | outward], dtype=np.longdouble))
    inner_edge = {
        "inward_count": int(np.sum(inward)),
        "inward_mass": float(np.sum(mass[inward], dtype=np.longdouble)),
        "outward_count": int(np.sum(outward)),
        "outward_mass": float(np.sum(mass[outward], dtype=np.longdouble)),
        "total_crossing_mass_over_initial_disk_mass": crossing_mass / disk["mass"],
        "definition": "Fixed-mass generator crossings of r=.5; not a face-integrated Eulerian flux",
    }
    log = (args.run / "log").read_text(errors="replace")
    energy_lines = [line for line in log.splitlines() if "[RK4E]" in line]
    if not energy_lines:
        raise ValueError("Missing RK4 energy ledger")
    last = energy_lines[-1]
    zero_injection = True
    for key in ("floor_cum", "de_cum", "sfl_cum"):
        match = re.search(rf"\b{key}=([0-9.eE+-]+)", last)
        zero_injection &= bool(match and float(match.group(1)) == 0.)
    energy_match = re.search(r"\bdEtot/\|Etot0\|=([0-9.eE+-]+)", last)
    energy_drift = float(energy_match.group(1)) if energy_match else float("nan")
    initial_volume_sum = float(np.sum(p0["volume"].astype(float), dtype=np.longdouble))
    endpoint_volume_sum = float(np.sum(p["volume"].astype(float), dtype=np.longdouble))
    checks = {"finite": bool(all(np.isfinite(q[field]).all() for q in (p0, p) for field in
                                 ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2"))),
              "initial_periodic_volume": abs(initial_volume_sum - 36) < 1e-6,
              "endpoint_periodic_volume": abs(endpoint_volume_sum - 36) < 1e-5,
              "mass_exact": bool(np.array_equal(p0["mass"], p["mass"])),
              "W0_exact": bool(np.all(p0["w2"] == 0) and np.all(p["w2"] == 0)),
              "positive_state": bool(np.all(p["volume"] > 0) and np.all(p["ie"] > 0)),
              "no_floor_injection": zero_injection,
              "total_energy_relative_below_1e_minus_3": bool(np.isfinite(energy_drift) and abs(energy_drift) < 1e-3),
              "disk_density_abs_below_0p005": disk["density_abs_mass_weighted"] < .005,
              "disk_entropy_abs_below_0p01": disk["entropy_abs_mass_weighted"] < .01,
              "disk_position_error_below_1e_minus_4": disk["position_error_mass_weighted"] < 1e-4,
              "disk_radial_speed_rms_below_0p01": disk["radial_speed_rms_over_circular"] < .01,
              "first_ring_density_abs_below_0p05": first_ring["density_abs_mass_weighted"] < .05,
              "first_ring_entropy_abs_below_0p01": first_ring["entropy_abs_mass_weighted"] < .01,
              "first_ring_radial_speed_rms_below_0p05": first_ring["radial_speed_rms_over_circular"] < .05,
              "inner_edge_crossing_mass_below_1e_minus_4_disk":
                  inner_edge["total_crossing_mass_over_initial_disk_mass"] < 1e-4}
    report = {"initial": str(args.initial), "endpoint": str(endpoint),
              "endpoint_sha256": digest(endpoint), "time": float(m["t"]),
              "n": n, "checks": checks, "initial_volume_sum": initial_volume_sum,
              "endpoint_volume_sum": endpoint_volume_sum,
              "mass_sum": float(np.sum(mass, dtype=np.longdouble)),
              "min_internal_energy": float(np.min(p["ie"])),
              "min_volume": float(np.min(p["volume"])),
              "max_density": float(np.max(p["den"])),
              "total_energy_relative_drift": energy_drift,
              "inner_disk_edge_cell_center_crossings": inner_edge,
              "bins": bins,
              "scope": "Corrected sparse fixed-mass W0 native CUDA+MPI+OpenMP early gate; does not establish long Kepler stability"}
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report, indent=2))
    if not all(checks.values()):
        raise SystemExit("Corrected sparse t=.825 physical gate failed")


if __name__ == "__main__":
    main()

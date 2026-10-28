"""Radial and tail diagnostics from an already accepted bounded Kepler state.

This is a read-only postprocessor. It does not alter a running simulation.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_native_ballistic_20261027 import disk_geometry
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
GM1 = float(np.float32(1.666667)) - 1.0
EDGES = (0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)


def weighted_mean(values, weights):
    return float(np.sum(values * weights, dtype=np.longdouble)
                 / np.sum(weights, dtype=np.longdouble))


def describe(mask, mass, entropy, density_residual, position_error,
             velocity_error, radial_speed):
    w = mass[mask]
    if not len(w):
        raise ValueError("Empty annulus")
    e = entropy[mask]
    d = density_residual[mask]
    return {
        "cells": int(mask.sum()),
        "mass": float(np.sum(w, dtype=np.longdouble)),
        "entropy_log_signed_mass_weighted": weighted_mean(e, w),
        "entropy_log_abs_mass_weighted": weighted_mean(np.abs(e), w),
        "entropy_log_abs_p99": float(np.percentile(np.abs(e), 99)),
        "entropy_log_abs_max": float(np.max(np.abs(e))),
        "density_native_minus_ballistic_abs_mass_weighted": weighted_mean(np.abs(d), w),
        "density_native_minus_ballistic_abs_p99": float(np.percentile(np.abs(d), 99)),
        "density_native_minus_ballistic_abs_max": float(np.max(np.abs(d))),
        "position_error_rms_mass_weighted": float(np.sqrt(weighted_mean(position_error[mask] ** 2, w))),
        "velocity_error_over_circular_rms_mass_weighted": float(np.sqrt(weighted_mean(velocity_error[mask] ** 2, w))),
        "radial_speed_over_circular_rms_mass_weighted": float(np.sqrt(weighted_mean(radial_speed[mask] ** 2, w))),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("accepted_receipt", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-failed-check", action="append", default=[],
                        help="Exploratory diagnostic only: require exactly this failed native gate")
    args = parser.parse_args()
    receipt = json.loads(args.accepted_receipt.read_text())
    failed = {key for key, passed in receipt["checks"].items() if not passed}
    if failed != set(args.allow_failed_check):
        raise ValueError(f"Native failed checks {sorted(failed)} differ from explicitly allowed {sorted(set(args.allow_failed_check))}")
    initial = Path(receipt["initial"])
    endpoint = Path(receipt["endpoint"])
    if digest(initial) != receipt["initial_sha256"] or digest(endpoint) != receipt["endpoint_sha256"]:
        raise ValueError("Native snapshot hash differs from accepted receipt")
    p0, _ = ordered(initial, N)
    p, _ = ordered(endpoint, N)
    x0 = p0["x"].astype(float) - 3.0
    y0 = p0["y"].astype(float) - 3.0
    r0 = np.hypot(x0, y0)
    disk = (r0 >= EDGES[0]) & (r0 < EDGES[-1])
    selected = np.flatnonzero(disk)
    time = float(receipt["time"])
    xy, ballistic_area, _ = disk_geometry(p0, time, selected)
    mass = p0["mass"][selected].astype(float)
    volume = p["volume"][selected].astype(float)
    volume0 = p0["volume"][selected].astype(float)
    entropy = (np.log(p["ie"][selected].astype(float) / p0["ie"][selected].astype(float))
               + GM1 * np.log(volume / volume0))
    density_residual = ballistic_area / volume - 1.0
    position_error = np.hypot(p["x"][selected].astype(float) - xy[selected, 0],
                              p["y"][selected].astype(float) - xy[selected, 1])
    radius = r0[selected]
    omega = (radius * radius + 0.01175 ** 2) ** (-0.75)
    angle = omega * time
    c, s = np.cos(angle), np.sin(angle)
    vx0 = p0["vx"][selected].astype(float)
    vy0 = p0["vy"][selected].astype(float)
    vx_expected = c * vx0 - s * vy0
    vy_expected = s * vx0 + c * vy0
    vx = p["vx"][selected].astype(float)
    vy = p["vy"][selected].astype(float)
    vk = radius * omega
    velocity_error = np.hypot(vx - vx_expected, vy - vy_expected) / vk
    x = p["x"][selected].astype(float) - 3.0
    y = p["y"][selected].astype(float) - 3.0
    radial_speed = np.abs(x * vx + y * vy) / np.maximum(np.hypot(x, y), 1e-30) / vk
    if not all(np.all(np.isfinite(a)) for a in
               (entropy, density_residual, position_error, velocity_error, radial_speed)):
        raise ValueError("Nonfinite diagnostic value")
    bands = []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        mask = (radius >= lo) & (radius < hi)
        bands.append({"r_min": lo, "r_max": hi,
                      **describe(mask, mass, entropy, density_residual,
                                 position_error, velocity_error, radial_speed)})
    whole = describe(np.ones(len(selected), dtype=bool), mass, entropy,
                     density_residual, position_error, velocity_error, radial_speed)
    if abs(whole["entropy_log_abs_mass_weighted"] - receipt["disk_entropy_log_abs_mass_weighted"]) > 1e-12:
        raise ValueError("Entropy metric does not reproduce parent audit")
    if abs(whole["density_native_minus_ballistic_abs_mass_weighted"]
           - receipt["disk_native_minus_ballistic_abs_mass_weighted"]) > 1e-12:
        raise ValueError("Density metric does not reproduce parent audit")
    result = {
        "accepted_receipt": str(args.accepted_receipt.resolve()),
        "initial_sha256": receipt["initial_sha256"],
        "endpoint_sha256": receipt["endpoint_sha256"],
        "time": time,
        "native_failed_checks": sorted(failed),
        "whole_disk": whole,
        "radial_bands": bands,
        "scope": "One-rank native bounded W0 state. Listed native failed checks remain failures; ballistic reference keeps exact circular site trajectories.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"time": time, "whole_disk": whole, "radial_bands": bands}, indent=2))


if __name__ == "__main__":
    main()

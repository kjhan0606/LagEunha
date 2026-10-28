#!/usr/bin/env python3
"""Recompute periodic-aware ballistic errors for the area-CFL t=2 checkpoint."""

import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot  # noqa: E402


ROOT = Path(
    "/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/"
    "area_cfl_apply_t1515_to65_a40_20261028"
)
INITIAL = ROOT / "khout.000000.dat"
ENDPOINT = ROOT / "khout.014793.dat"
OUTPUT = Path(__file__).with_name(
    "KEPLER_Q1_AREA_CFL_T2_PERIODIC_RADIAL_AUDIT_2026-10-28.json"
)
N = 384 * 384
PERIOD = 6.0
CENTER = 3.0
EPS = 0.01175
MASK = np.uint64((1 << 56) - 1)


def load(path, expected_time):
    data, meta = read_snapshot(path)
    if len(data) != N or abs(float(meta["t"]) - expected_time) > 1e-10:
        raise ValueError(f"Unexpected snapshot shape/time: {path}, {meta['t']}")
    ids = data["u4if"].astype(np.uint64) & MASK
    order = np.argsort(ids)
    if not np.array_equal(ids[order], np.arange(N)):
        raise ValueError(f"Incomplete or duplicated IDs: {path}")
    return np.asarray(data[order]), float(meta["t"])


def weighted_mean(values, weights):
    return float(np.sum(values * weights, dtype=np.longdouble)
                 / np.sum(weights, dtype=np.longdouble))


def main():
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUTPUT}")
    initial, t0 = load(INITIAL, 0.0)
    final, time = load(ENDPOINT, 2.0001076440463654)
    if not np.array_equal(initial["mass"], final["mass"]):
        raise ValueError("Mass is not bitwise fixed by ID")

    x0 = initial["x"].astype(float) - CENTER
    y0 = initial["y"].astype(float) - CENTER
    r0 = np.hypot(x0, y0)
    omega = (r0 * r0 + EPS * EPS) ** (-0.75)
    phase = time * omega
    ballistic_x = CENTER + x0 * np.cos(phase) - y0 * np.sin(phase)
    ballistic_y = CENTER + x0 * np.sin(phase) + y0 * np.cos(phase)

    # Compare on the periodic 6x6 torus so a boundary crossing is not counted
    # as an O(box-size) physical position error.
    dx = (final["x"].astype(float) - ballistic_x + PERIOD / 2) % PERIOD - PERIOD / 2
    dy = (final["y"].astype(float) - ballistic_y + PERIOD / 2) % PERIOD - PERIOD / 2
    error = np.hypot(dx, dy)
    mass = initial["mass"].astype(float)
    bins = []
    for lo, hi in ((0.0, 0.1), (0.1, 0.5), (0.5, 0.75), (0.75, 1.0),
                   (1.0, 2.0), (0.5, 2.0), (2.0, 2.35), (2.35, 2.75),
                   (2.75, 5.0)):
        select = (r0 >= lo) & (r0 < hi)
        bins.append({
            "initial_radius": [lo, hi],
            "cells": int(select.sum()),
            "mass_weighted_periodic_ballistic_position_error": weighted_mean(
                error[select], mass[select]),
        })
    result = {
        "initial": str(INITIAL),
        "endpoint": str(ENDPOINT),
        "initial_time": t0,
        "endpoint_time": time,
        "domain_period": PERIOD,
        "fixed_mass_bitwise_by_id": True,
        "method": "initial softened-Kepler ballistic orbit; wrapped Cartesian displacement on 6x6 periodic torus",
        "bins": bins,
        "note": "The unwrapped outer-bin position error in the companion long-gate receipt is not a physical periodic displacement.",
    }
    with OUTPUT.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))
    print(f"output={OUTPUT.resolve()}")


if __name__ == "__main__":
    main()

"""Compare per-cell orbital angular momentum in a same-IC Kepler pair."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


BINS = ((.5, .52), (.52, .55), (.55, .65), (.5, .65),
        (.65, .75), (.75, 1.), (1., 2.), (.5, 2.))


def weighted_mean_abs(values, masses):
    return float(np.sum(np.abs(values) * masses, dtype=np.longdouble)
                 / np.sum(masses, dtype=np.longdouble))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("w0_receipt", type=Path)
    parser.add_argument("weighted_receipt", type=Path)
    parser.add_argument("--labels", nargs=2, default=("w0", "weighted"))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    first, second = args.labels
    if first == second:
        raise ValueError("Labels must differ")
    receipts = {name: json.loads(path.read_text()) for name, path in
                ((first, args.w0_receipt), (second, args.weighted_receipt))}
    if receipts[first]["initial_sha256"] != receipts[second]["initial_sha256"]:
        raise ValueError("Initial snapshots differ")
    if receipts[first]["time"] != receipts[second]["time"]:
        raise ValueError("Endpoint times differ")
    if receipts[first]["n"] != receipts[second]["n"]:
        raise ValueError("Site counts differ")

    n = receipts[first]["n"]
    initial = Path(receipts[first]["initial"])
    if digest(initial) != receipts[first]["initial_sha256"]:
        raise ValueError("Initial snapshot digest mismatch")
    p0, meta0 = ordered(initial, n)
    if abs(meta0["t"]) > 1e-12:
        raise ValueError("Initial snapshot is not at t=0")
    x0, y0 = p0["x"].astype(float) - 3., p0["y"].astype(float) - 3.
    r0 = np.hypot(x0, y0)
    mass = p0["mass"].astype(float)
    ell0 = x0 * p0["vy"].astype(float) - y0 * p0["vx"].astype(float)
    if np.any(ell0[(r0 >= .5) & (r0 < 2.)] <= 0):
        raise ValueError("Expected prograde initial disk")

    report = {"initial_sha256": receipts[first]["initial_sha256"],
              "time": receipts[first]["time"], "n": n, "bins": {}}
    for name, receipt in receipts.items():
        endpoint = Path(receipt["endpoint"])
        if digest(endpoint) != receipt["endpoint_sha256"]:
            raise ValueError(f"{name} endpoint digest mismatch")
        p, meta = ordered(endpoint, n)
        if abs(meta["t"] - report["time"]) > 1e-12:
            raise ValueError(f"{name} endpoint time mismatch")
        if not np.array_equal(p["mass"], p0["mass"]):
            raise ValueError(f"{name} cell masses changed")
        x, y = p["x"].astype(float) - 3., p["y"].astype(float) - 3.
        ell = x * p["vy"].astype(float) - y * p["vx"].astype(float)
        rel = ell / ell0 - 1.
        rows = {}
        for lo, hi in BINS:
            take = (r0 >= lo) & (r0 < hi)
            w = mass[take]
            delta = rel[take]
            l0 = np.sum(w * ell0[take], dtype=np.longdouble)
            l1 = np.sum(w * ell[take], dtype=np.longdouble)
            rows[f"{lo:g}:{hi:g}"] = {
                "cells": int(take.sum()),
                "mean_abs_relative_per_cell": weighted_mean_abs(delta, w),
                "p99_abs_relative_per_cell": float(np.quantile(np.abs(delta), .99)),
                "max_abs_relative_per_cell": float(np.max(np.abs(delta))),
                "signed_total_relative": float(l1 / l0 - 1.),
            }
        report["bins"][name] = rows
    report["scope"] = ("Initial-ID bins; orbital centre-of-mass angular momentum "
                       "only. Global signed cancellation can hide local torque. "
                       "No cell spin angular momentum or material-face flux measured.")
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    for lo, hi in BINS:
        key = f"{lo:g}:{hi:g}"
        left = report["bins"][first][key]
        right = report["bins"][second][key]
        print(key, "mean_abs_relative", left["mean_abs_relative_per_cell"],
              right["mean_abs_relative_per_cell"], "ratio",
              left["mean_abs_relative_per_cell"] / right["mean_abs_relative_per_cell"])


if __name__ == "__main__":
    main()

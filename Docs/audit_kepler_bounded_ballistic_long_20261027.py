"""Exact-orbit Voronoi geometry of both native bounded Kepler initial ICs."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_native_ballistic_20261027 import disk_geometry, weighted_abs
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N = 55351
EPS = .01175
T_NOMINAL = 1000 * 2 * np.pi * .5 ** 1.5
T_SOFTENED = 1000 * 2 * np.pi * (.5 ** 2 + EPS ** 2) ** .75
TIMES = (6.5, 600., T_NOMINAL, T_SOFTENED)
BANDS = ((.5, .75), (.75, 1.), (1., 2.), (.5, 2.))


def one_case(label, receipt_path):
    receipt = json.loads(receipt_path.read_text())
    if not all(receipt["checks"].values()):
        raise ValueError(f"{label}: native parent audit failed")
    initial = Path(receipt["initial"])
    if digest(initial) != receipt["initial_sha256"] or abs(snapshot_time(initial, N)) > 1e-12:
        raise ValueError(f"{label}: native initial snapshot changed")
    p0, _ = ordered(initial, N)
    if not np.all(p0["w2"] == 0):
        raise ValueError(f"{label}: W must be zero")
    r0 = np.hypot(p0["x"].astype(float) - 3., p0["y"].astype(float) - 3.)
    selected = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    if len(selected) != 49038:
        raise ValueError(f"{label}: disk site count changed")
    mass = p0["mass"][selected].astype(float)
    volume0 = p0["volume"][selected].astype(float)
    times = []
    for t in TIMES:
        _, area, safety = disk_geometry(p0, float(t), selected)
        density_change = volume0 / area - 1.
        bands = []
        for lo, hi in BANDS:
            mask = (r0[selected] >= lo) & (r0[selected] < hi)
            err = np.abs(density_change[mask])
            bands.append({
                "r_min": lo, "r_max": hi, "cells": int(np.sum(mask)),
                "density_abs_mass_weighted": weighted_abs(density_change[mask], mass[mask]),
                "density_abs_p99": float(np.percentile(err, 99)),
                "density_abs_max": float(np.max(err)),
            })
        times.append({"time": float(t), "minimum_periodic_image_safety_ratio": safety,
                      "bands": bands})
    return {"label": label, "receipt": str(receipt_path.resolve()),
            "initial": str(initial.resolve()), "initial_sha256": receipt["initial_sha256"],
            "disk_sites": int(len(selected)), "times": times}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("aligned_accepted_receipt", type=Path)
    parser.add_argument("jitter_accepted_receipt", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    aligned = one_case("aligned", args.aligned_accepted_receipt)
    jitter = one_case("5% independent site jitter", args.jitter_accepted_receipt)
    ratios = []
    for a, j in zip(aligned["times"], jitter["times"]):
        if a["time"] != j["time"]:
            raise ValueError("Mismatched orbit-only time")
        ratios.append({
            "time": a["time"],
            "jitter_over_aligned_disk_density_abs":
                j["bands"][-1]["density_abs_mass_weighted"]
                / a["bands"][-1]["density_abs_mass_weighted"],
        })
    result = {
        "nominal_unsoftened_1000_inner_orbit_time": float(T_NOMINAL),
        "actual_softened_1000_inner_orbit_time": float(T_SOFTENED),
        "inner_reference_radius": .5,
        "epsilon": EPS,
        "cases": [aligned, jitter],
        "ratios": ratios,
        "scope": "Exact circular point orbits and ordinary W=0 Voronoi areas for the two accepted bounded ICs; fixed initial masses. No native hydro, wall response, entropy, or proof of long-time stability.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"times": list(TIMES),
                      "aligned_disk_abs": [t["bands"][-1]["density_abs_mass_weighted"] for t in aligned["times"]],
                      "jitter_disk_abs": [t["bands"][-1]["density_abs_mass_weighted"] for t in jitter["times"]],
                      "ratios": ratios}, indent=2))


if __name__ == "__main__":
    main()

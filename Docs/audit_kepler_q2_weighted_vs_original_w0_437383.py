"""Matched native weighted step 437383 versus the original W0 CUDA binary."""

import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


BASE = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
ORIGINAL = BASE / "faceforce_pair_t6501_a10_v1"
WEIGHTED = BASE / "gpu_energy_neutral_weighted_t6501_guard_a10_v2"
OUTPUT = Path(__file__).with_name("KEPLER_Q2_WEIGHTED_VS_ORIGINAL_W0_437383.json")
N = 55351
GAMMA = float(np.float32(1.666667))


def read_receipt(root, label):
    d = json.loads((root / label / "audit_t6501.json").read_text())
    if (not all(d["checks"].values()) or d["start_time"] != 6.5
            or d["endpoint_time"] != 6.501):
        raise ValueError(f"{label}: unaccepted endpoint in {root}")
    for key, time in (("restart", 6.5), ("endpoint", 6.501)):
        path = Path(d[key])
        if (path.parent != root / label or digest(path) != d[key + "_sha256"]
                or abs(snapshot_time(path, N)-time)>1e-10):
            raise ValueError(f"{label}: changed {key} in {root}")
    return d


def metrics(start, end, chosen):
    mass = start["mass"][chosen].astype(float)
    volume = end["volume"][chosen].astype(float)/start["volume"][chosen].astype(float)-1.
    density = np.log(end["den"][chosen].astype(float)/start["den"][chosen].astype(float))
    pressure = np.log(end["pressure"][chosen].astype(float)/start["pressure"][chosen].astype(float))
    entropy = pressure-GAMMA*density
    end_pressure = end["pressure"][chosen].astype(float)
    mean_pressure = np.average(end_pressure, weights=mass)
    return {
        "cells": int(np.count_nonzero(chosen)),
        "volume_change_abs_mass_mean": float(np.average(np.abs(volume), weights=mass)),
        "density_log_change_abs_mass_mean": float(np.average(np.abs(density), weights=mass)),
        "pressure_log_change_rms": float(np.sqrt(np.average(pressure**2, weights=mass))),
        "pressure_std_over_mean": float(np.sqrt(np.average(
            (end_pressure-mean_pressure)**2, weights=mass))/mean_pressure),
        "entropy_log_change_abs_mass_mean": float(np.average(np.abs(entropy), weights=mass)),
    }


def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    results = {
        "time_start": 6.5, "time_end": 6.501,
        "original_W0_run": str(ORIGINAL), "weighted_run": str(WEIGHTED),
        "original_W0_binary_sha256": "3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e",
        "weighted_binary_sha256": "e90a287d25b0d82705f70ad645b5eadade7a18cd3aebff1622606960a333e7a1",
        "cases": {},
        "scope": "Matched one-step original CUDA W0 versus experimental weighted native RK4. "
                 "Local moving-face first law and long-time Kepler stability remain open.",
    }
    for label in ("aligned", "jitter"):
        z = read_receipt(ORIGINAL, label)
        w = read_receipt(WEIGHTED, label)
        if z["restart_sha256"] != w["restart_sha256"]:
            raise ValueError(f"{label}: different physical parent")
        start, _ = ordered(Path(z["restart"]), N)
        original, _ = ordered(Path(z["endpoint"]), N)
        weighted, _ = ordered(Path(w["endpoint"]), N)
        if (not np.array_equal(start["mass"], original["mass"])
                or not np.array_equal(start["mass"], weighted["mass"])
                or not np.all(original["w2"] == 0)):
            raise ValueError(f"{label}: fixed-mass W0 reduction failed")
        radius = np.hypot(start["x"].astype(float)-3., start["y"].astype(float)-3.)
        bands = {}
        for lo, hi in ((.5, .75), (.75, 1.), (1., 2.), (.5, 2.)):
            chosen = (radius>=lo)&(radius<hi)
            a, b = metrics(start, original, chosen), metrics(start, weighted, chosen)
            bands[f"{lo:g}-{hi:g}"] = {
                "original_W0": a, "weighted": b,
                "weighted_over_W0_volume_change":
                    b["volume_change_abs_mass_mean"]/a["volume_change_abs_mass_mean"],
            }
        results["cases"][label] = {
            "parent_sha256": z["restart_sha256"],
            "original_W0_endpoint_sha256": z["endpoint_sha256"],
            "weighted_endpoint_sha256": w["endpoint_sha256"],
            "original_W0_energy_relative_change": z["energy_relative_change_from_restart"],
            "weighted_energy_relative_change": w["energy_relative_change_from_restart"],
            "weighted_max_weight": w["weight_max"],
            "bands": bands,
        }
    with OUTPUT.open("x") as output:
        json.dump(results, output, indent=2, allow_nan=False)
        output.write("\n")
    print(json.dumps({label: {
        "whole_disk_volume_ratio": case["bands"]["0.5-2"]["weighted_over_W0_volume_change"],
        "inner_volume_ratio": case["bands"]["0.5-0.75"]["weighted_over_W0_volume_change"],
    } for label, case in results["cases"].items()}, indent=2))


if __name__ == "__main__":
    main()

"""Compare native Kepler pressure with adiabatic pressure from ballistic cells."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_native_ballistic_20261027 import disk_geometry
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
GAMMA = float(np.float32(1.666667))


def weighted_mean(values, weights):
    return float(np.sum(values * weights, dtype=np.longdouble)
                 / np.sum(weights, dtype=np.longdouble))


def describe(mask, mass, actual, predicted):
    a, p, w = actual[mask], predicted[mask], mass[mask]
    ma, mp = weighted_mean(a, w), weighted_mean(p, w)
    da, dp = a - ma, p - mp
    covariance = weighted_mean(da * dp, w)
    variance_a = weighted_mean(da * da, w)
    variance_p = weighted_mean(dp * dp, w)
    return {
        "cells": int(np.count_nonzero(mask)),
        "native_pressure_relative_abs_mass_weighted": weighted_mean(np.abs(a), w),
        "ballistic_adiabatic_pressure_relative_abs_mass_weighted": weighted_mean(np.abs(p), w),
        "native_minus_ballistic_pressure_relative_abs_mass_weighted": weighted_mean(np.abs(a-p), w),
        "native_pressure_relative_signed_mass_weighted": ma,
        "ballistic_pressure_relative_signed_mass_weighted": mp,
        "pressure_relative_pearson_mass_weighted": float(covariance / np.sqrt(variance_a*variance_p))
        if variance_a > 0 and variance_p > 0 else None,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accepted_receipt", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    receipt = json.loads(args.accepted_receipt.read_text())
    if len(receipt["checks"]) != 17 or not all(receipt["checks"].values()):
        raise ValueError("Parent native state failed a 17-check gate")
    initial, endpoint = Path(receipt["initial"]), Path(receipt["endpoint"])
    if digest(initial) != receipt["initial_sha256"] or digest(endpoint) != receipt["endpoint_sha256"]:
        raise ValueError("Native snapshot changed from accepted receipt")
    p0, meta0 = ordered(initial, N)
    p, meta = ordered(endpoint, N)
    time = float(receipt["time"])
    if abs(meta0["t"]) > 1e-12 or abs(meta["t"]-time) > 1e-10:
        raise ValueError("Native checkpoint time mismatch")
    r0 = np.hypot(p0["x"].astype(float)-3., p0["y"].astype(float)-3.)
    selected = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    _, area_ballistic, _ = disk_geometry(p0, time, selected)
    mass = p0["mass"][selected].astype(float)
    volume0 = p0["volume"][selected].astype(float)
    pressure0 = p0["pressure"][selected].astype(float)
    pressure = p["pressure"][selected].astype(float)
    actual = pressure/pressure0-1.
    predicted = (volume0/area_ballistic)**GAMMA-1.
    if not all(np.all(np.isfinite(q)) for q in (actual, predicted)):
        raise ValueError("Nonfinite pressure diagnostic")
    radius = r0[selected]
    result = {
        "accepted_receipt": str(args.accepted_receipt.resolve()),
        "initial_sha256": receipt["initial_sha256"],
        "endpoint_sha256": receipt["endpoint_sha256"],
        "time": time,
        "gamma": GAMMA,
        "whole_disk": describe(np.ones(len(selected), dtype=bool), mass, actual, predicted),
        "inner_r05_r075": describe((radius >= .5) & (radius < .75), mass, actual, predicted),
        "scope": ("Adiabatic pressure predicted from exact circular-site Voronoi area and each site's initial "
                  "pressure. Agreement measures geometric pressure excitation, not the native per-face pressure force."),
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

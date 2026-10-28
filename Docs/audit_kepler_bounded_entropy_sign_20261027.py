"""Resolve signed Kepler entropy error into positive and negative parts.

This reads an accepted native snapshot pair. It never modifies a simulation.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
GM1 = float(np.float32(1.666667)) - 1.0
EDGES = (0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
REQUIRED_CHECKS = {
    "finite_positive", "fixed_mass_exact", "W0_exact",
    "initial_area_tiles_polygon", "endpoint_area_tiles_polygon",
    "density_and_EOS", "zero_injection", "total_energy", "disk_density",
    "disk_entropy", "disk_orbit", "disk_ballistic_residual",
    "disk_radial_velocity", "disk_angular_momentum", "inside_polygon_incircle",
    "no_wall_crossing", "fixed_polygon_enabled",
}


def describe(entropy, mass):
    total = np.sum(mass, dtype=np.longdouble)
    if not np.isfinite(total) or total <= 0:
        raise ValueError("Empty or invalid annulus")

    def mean(value):
        return float(np.sum(value * mass, dtype=np.longdouble) / total)

    positive = np.maximum(entropy, 0.0)
    negative = np.maximum(-entropy, 0.0)
    result = {
        "cells": len(entropy),
        "mass": float(total),
        "entropy_log_signed_mass_weighted": mean(entropy),
        "entropy_log_abs_mass_weighted": mean(np.abs(entropy)),
        "entropy_log_positive_component_mass_weighted": mean(positive),
        "entropy_log_negative_component_mass_weighted": mean(negative),
        "positive_cell_fraction": float(np.count_nonzero(entropy > 0) / len(entropy)),
        "positive_mass_fraction": mean((entropy > 0).astype(float)),
        "negative_cell_fraction": float(np.count_nonzero(entropy < 0) / len(entropy)),
        "negative_mass_fraction": mean((entropy < 0).astype(float)),
        "entropy_log_p01": float(np.quantile(entropy, 0.01)),
        "entropy_log_median": float(np.median(entropy)),
        "entropy_log_p99": float(np.quantile(entropy, 0.99)),
    }
    if abs(result["entropy_log_signed_mass_weighted"] -
           (result["entropy_log_positive_component_mass_weighted"] -
            result["entropy_log_negative_component_mass_weighted"])) > 1e-15:
        raise ValueError("Signed entropy decomposition failed")
    if abs(result["entropy_log_abs_mass_weighted"] -
           (result["entropy_log_positive_component_mass_weighted"] +
            result["entropy_log_negative_component_mass_weighted"])) > 1e-15:
        raise ValueError("Absolute entropy decomposition failed")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("accepted_receipt", type=Path)
    parser.add_argument("--allow-failed-check", action="append", default=[],
                        help="Diagnose the signed entropy state despite these explicitly named failed gates")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads(args.accepted_receipt.read_text())
    failed = {key for key, value in receipt["checks"].items() if not value}
    if (not REQUIRED_CHECKS.issubset(receipt["checks"]) or
            failed != set(args.allow_failed_check)):
        raise ValueError("Parent native receipt has an unlisted failed gate or is missing a required gate")
    initial, endpoint = Path(receipt["initial"]), Path(receipt["endpoint"])
    if digest(initial) != receipt["initial_sha256"] or digest(endpoint) != receipt["endpoint_sha256"]:
        raise ValueError("Native checkpoint hash differs from accepted receipt")
    p0, m0 = ordered(initial, N)
    p, m = ordered(endpoint, N)
    if abs(m0["t"]) > 1e-12 or abs(m["t"] - receipt["time"]) > 1e-10:
        raise ValueError("Native checkpoint time differs from accepted receipt")
    r0 = np.hypot(p0["x"].astype(float) - 3.0, p0["y"].astype(float) - 3.0)
    selected = (r0 >= EDGES[0]) & (r0 < EDGES[-1])
    entropy = (np.log(p["ie"][selected].astype(float) / p0["ie"][selected].astype(float))
               + GM1 * np.log(p["volume"][selected].astype(float) / p0["volume"][selected].astype(float)))
    mass = p0["mass"][selected].astype(float)
    radius = r0[selected]
    if not np.all(np.isfinite(entropy)):
        raise ValueError("Nonfinite entropy diagnostic")
    whole = describe(entropy, mass)
    if abs(whole["entropy_log_abs_mass_weighted"] - receipt["disk_entropy_log_abs_mass_weighted"]) > 1e-12:
        raise ValueError("Entropy diagnostic differs from native receipt")
    bands = []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        mask = (radius >= lo) & (radius < hi)
        bands.append({"r_min": lo, "r_max": hi, **describe(entropy[mask], mass[mask])})
    result = {
        "accepted_receipt": str(args.accepted_receipt.resolve()),
        "initial_sha256": receipt["initial_sha256"],
        "endpoint_sha256": receipt["endpoint_sha256"],
        "time": receipt["time"],
        "native_failed_checks": sorted(failed),
        "whole_disk": whole,
        "radial_bands": bands,
        "scope": ("Signed log-entropy relative to each site's initial state; positive and negative "
                  "components are diagnostic errors, not separately measured irreversible heating."),
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"time": result["time"], "whole_disk": whole,
                      "inner_disk": bands[0]}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

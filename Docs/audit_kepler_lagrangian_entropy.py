#!/usr/bin/env python3
"""Audit fixed-ID Kepler disk and ambient cohorts across native snapshots."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot


MASK56 = np.uint64((1 << 56) - 1)
BANDS = (
    ("central", 0.0, 0.5),
    ("disk_inner", 0.5, 1.0),
    ("disk_outer", 1.0, 2.0),
    ("disk", 0.5, 2.0),
    ("edge_transition", 2.0, 2.5),
    ("ambient", 2.5, 4.3),
)


def ordered_snapshot(path: Path) -> tuple[np.ndarray, dict, np.ndarray]:
    particles, meta = read_snapshot(path)
    ids = particles["u4if"].astype(np.uint64) & MASK56
    order = np.argsort(ids)
    return particles[order], meta, ids[order]


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.dot(weights, values) / np.sum(weights))


def audit(initial_path: Path, target_path: Path, gamma: float,
          center: tuple[float, float]) -> dict:
    initial, initial_meta, initial_ids = ordered_snapshot(initial_path)
    target, target_meta, target_ids = ordered_snapshot(target_path)
    if not np.array_equal(initial_ids, target_ids):
        raise ValueError(f"Cell ID sets differ: {initial_path} vs {target_path}")

    mass0 = np.asarray(initial["mass"], dtype=np.float64)
    mass = np.asarray(target["mass"], dtype=np.float64)
    if not np.array_equal(mass0, mass):
        raise ValueError(f"Fixed-mass invariant failed: {initial_path} vs {target_path}")

    rho0 = np.asarray(initial["den"], dtype=np.float64)
    rho = np.asarray(target["den"], dtype=np.float64)
    pressure0 = np.asarray(initial["pressure"], dtype=np.float64)
    pressure = np.asarray(target["pressure"], dtype=np.float64)
    ie0 = np.asarray(initial["ie"], dtype=np.float64)
    ie = np.asarray(target["ie"], dtype=np.float64)
    if min(np.min(rho0), np.min(rho), np.min(pressure0), np.min(pressure)) <= 0:
        raise ValueError("Entropy proxy requires positive density and pressure")

    x0 = np.asarray(initial["x"], dtype=np.float64) - center[0]
    y0 = np.asarray(initial["y"], dtype=np.float64) - center[1]
    x = np.asarray(target["x"], dtype=np.float64) - center[0]
    y = np.asarray(target["y"], dtype=np.float64) - center[1]
    r0 = np.hypot(x0, y0)
    r = np.hypot(x, y)
    delta_ln_rho = np.log(rho / rho0)
    delta_ln_pressure = np.log(pressure / pressure0)
    delta_ln_u = np.log((ie / mass) / (ie0 / mass0))
    delta_ln_entropy = delta_ln_pressure - gamma * delta_ln_rho
    mass_sum = float(mass.sum())

    bands = {}
    for name, lower, upper in BANDS:
        selected = (r0 >= lower) & (r0 < upper)
        if not np.any(selected):
            bands[name] = {
                "initial_radius": [lower, upper],
                "cell_count": 0,
                "empty": True,
            }
            continue
        weights = mass[selected]
        weights /= weights.sum()
        radial_valid = r0[selected] > 1.0e-8
        radial_drift = ((r[selected][radial_valid] - r0[selected][radial_valid]) /
                        r0[selected][radial_valid])
        bands[name] = {
            "initial_radius": [lower, upper],
            "cell_count": int(np.count_nonzero(selected)),
            "mass_fraction": float(mass[selected].sum() / mass_sum),
            "internal_energy_sum_ratio": float(ie[selected].sum() / ie0[selected].sum()),
            "mass_weighted_mean_dln_density": weighted_mean(delta_ln_rho[selected], weights),
            "mass_weighted_mean_dln_pressure": weighted_mean(delta_ln_pressure[selected], weights),
            "mass_weighted_mean_dln_specific_internal_energy": weighted_mean(delta_ln_u[selected], weights),
            "mass_weighted_mean_dln_entropy_proxy": weighted_mean(delta_ln_entropy[selected], weights),
            "mass_weighted_mean_abs_dln_entropy_proxy": weighted_mean(np.abs(delta_ln_entropy[selected]), weights),
            "dln_pressure_minus_gamma_dln_density": weighted_mean(
                delta_ln_pressure[selected] - gamma * delta_ln_rho[selected], weights),
            "dln_u_minus_gamma_minus_1_dln_density": weighted_mean(
                delta_ln_u[selected] - (gamma - 1.0) * delta_ln_rho[selected], weights),
            "rms_fractional_radius_drift": float(np.sqrt(np.mean(radial_drift ** 2))),
            "p99_abs_fractional_radius_drift": float(np.quantile(np.abs(radial_drift), 0.99)),
            "radius_drift_excluded_near_center_count": int(
                np.count_nonzero(~radial_valid)),
            "fraction_still_in_initial_radial_band": float(np.mean(
                (r[selected] >= lower) & (r[selected] < upper))),
        }

    return {
        "initial_snapshot": str(initial_path.resolve()),
        "target_snapshot": str(target_path.resolve()),
        "initial_time": float(initial_meta["t"]),
        "target_time": float(target_meta["t"]),
        "cell_count": int(len(mass)),
        "gamma": gamma,
        "center": list(center),
        "fixed_mass_bitwise_equal_by_cell_id": True,
        "total_mass": mass_sum,
        "global_internal_energy_sum_ratio": float(ie.sum() / ie0.sum()),
        "initial_radius_cohorts": bands,
        "interpretation": (
            "Cohorts are selected by initial radius and tracked by cell ID. "
            "The entropy proxy is ln(P/rho^gamma); the two residuals compare "
            "changes with ideal-gas adiabatic compression."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("targets", nargs="+", type=Path)
    parser.add_argument("--gamma", type=float, default=5.0 / 3.0)
    parser.add_argument("--center", nargs=2, type=float, default=(3.0, 3.0))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = {
        "method": "fixed-mass Lagrangian radial-cohort audit",
        "results": [audit(args.initial, target, args.gamma, tuple(args.center))
                    for target in args.targets],
    }
    serialized = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as stream:
            stream.write(serialized)
    print(serialized, end="")


if __name__ == "__main__":
    main()

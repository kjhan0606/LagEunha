"""Split saved coarse/fine relative-pressure differences into rho and u parts.

Read-only analysis of already accepted t=.05 checkpoints; no provider, RHS,
or new time step. The split is an algebraic identity, not causal attribution.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from compare_full_moment_checkpoint_errors import load


def weighted_rms(values, weights):
    values = np.asarray(values, dtype=np.longdouble)
    weights = np.asarray(weights, dtype=np.longdouble)
    return float(np.sqrt(np.sum(weights * values * values) / np.sum(weights)))


def weighted_corr(a, b, weights):
    a = np.asarray(a, dtype=np.longdouble)
    b = np.asarray(b, dtype=np.longdouble)
    weights = np.asarray(weights, dtype=np.longdouble)
    weights = weights / weights.sum()
    am = np.sum(weights * a)
    bm = np.sum(weights * b)
    av = np.sum(weights * (a - am) ** 2)
    bv = np.sum(weights * (b - bm) ** 2)
    if av == 0 or bv == 0:
        return None
    return float(np.sum(weights * (a - am) * (b - bm)) / np.sqrt(av * bv))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=Path(
        "/gpfs/kjhan/LagForce/Codex_full_disk_continuous_two_dt_to_t1_20261019"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"refusing to overwrite {args.output}")

    coarse, coarse_report, coarse_source = load(
        args.run_root / "coarse/checkpoints/coarse_0010/physics")
    fine, fine_report, fine_source = load(
        args.run_root / "fine/checkpoints/fine_0020/physics")
    if coarse_report["time"] != fine_report["time"] or coarse_report["time"] != 0.05:
        raise SystemExit("expected matching t=.05 checkpoints")
    if not np.array_equal(coarse["ids"], fine["ids"]):
        raise SystemExit("cell IDs differ")
    if not np.array_equal(coarse["actual_mass"], fine["actual_mass"]):
        raise SystemExit("fixed cell masses differ")

    # Ratios use each checkpoint's own steady reference on its actual cell.
    rho_c = (coarse["actual_density"] / coarse["reference_density"]).astype(np.longdouble)
    rho_f = (fine["actual_density"] / fine["reference_density"]).astype(np.longdouble)
    u_c = (coarse["actual_specific_thermal"]
           / coarse["reference_specific_thermal"]).astype(np.longdouble)
    u_f = (fine["actual_specific_thermal"]
           / fine["reference_specific_thermal"]).astype(np.longdouble)
    p_c = (coarse["actual_pressure"] / coarse["reference_pressure"]).astype(np.longdouble)
    p_f = (fine["actual_pressure"] / fine["reference_pressure"]).astype(np.longdouble)

    pressure_delta = p_f - p_c
    density_term = (rho_f - rho_c) * (u_f + u_c) / 2
    thermal_term = (u_f - u_c) * (rho_f + rho_c) / 2
    identity_residual = pressure_delta - density_term - thermal_term
    mass = coarse["actual_mass"].astype(np.longdouble)
    regions = {}
    for name in ("all", "cavity", "inner_contact", "disk_core"):
        if not np.array_equal(coarse[name], fine[name]):
            raise SystemExit(f"region mask differs: {name}")
        mask = coarse[name]
        w = mass[mask]
        d = pressure_delta[mask]
        x = density_term[mask]
        y = thermal_term[mask]
        regions[name] = {
            "cells": int(mask.sum()),
            "pressure_error_difference_RMS": weighted_rms(d, w),
            "density_contribution_RMS": weighted_rms(x, w),
            "specific_thermal_contribution_RMS": weighted_rms(y, w),
            "density_thermal_contribution_correlation": weighted_corr(x, y, w),
            "identity_residual_RMS": weighted_rms(identity_residual[mask], w),
            "identity_residual_maximum": float(np.max(np.abs(identity_residual[mask]))),
        }

    result = {
        "status": "completed_saved_pressure_error_product_decomposition",
        "time": 0.05,
        "equation": (
            "delta(P/P_ref) = delta(rho/rho_ref)*mean(u/u_ref) "
            "+ delta(u/u_ref)*mean(rho/rho_ref); each branch uses its own "
            "moving-cell steady reference"
        ),
        "regions": regions,
        "source_coarse": coarse_source,
        "source_fine": fine_source,
        "physical_providers": 0,
        "accepted_time_steps": 0,
        "limits": (
            "Exact algebraic decomposition of the difference between two "
            "checkpoint-relative pressure errors. It does not identify causal "
            "force or heating sources and does not establish convergence order."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "status": result["status"],
        "time": result["time"],
        "regions": {
            name: {
                key: row[key]
                for key in (
                    "pressure_error_difference_RMS",
                    "density_contribution_RMS",
                    "specific_thermal_contribution_RMS",
                    "density_thermal_contribution_correlation",
                    "identity_residual_maximum",
                )
            }
            for name, row in regions.items()
        },
    }, sort_keys=True))


if __name__ == "__main__":
    main()

"""Locate coarse/fine t=.05 temporal error contributions in saved checkpoints.

Read-only: no provider, RHS, time step, or checkpoint mutation. Errors are the
existing physical-error differences, mass weighted on the identical cell IDs.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from compare_full_moment_checkpoint_errors import load


def concentration(ids, radius_a, radius_b, error_a, error_b, delta, mass, mask):
    weighted = mass[mask].astype(np.longdouble)
    values = delta[mask].astype(np.longdouble)
    contribution = weighted * values * values
    total = contribution.sum()
    order_local = np.argsort(contribution)[::-1]
    order = np.flatnonzero(mask)[order_local]
    cumulative = np.cumsum(contribution[order_local]) / total
    radius_mid = 0.5 * (radius_a[mask] + radius_b[mask])
    radial_bands = {}
    for width in (0.005, 0.01, 0.02, 0.04):
        band = np.abs(radius_mid - 0.5) < width
        radial_bands[f"abs_r_minus_0.5_lt_{width:g}"] = {
            "cells": int(band.sum()),
            "contribution_fraction": float(contribution[band].sum() / total),
        }

    def count_for(fraction):
        return int(np.searchsorted(cumulative, fraction, side="left") + 1)

    top = []
    for rank, index in enumerate(order[:10]):
        top.append({
            "id": int(ids[index]),
            "radius_coarse": float(radius_a[index]),
            "radius_fine": float(radius_b[index]),
            "error_coarse": float(error_a[index]),
            "error_fine": float(error_b[index]),
            "delta": float(delta[index]),
            "contribution_fraction": float(contribution[order_local[rank]] / total),
        })
    return {
        "cells": int(mask.sum()),
        "weighted_RMS": float(np.sqrt(total / weighted.sum())),
        "cells_for_50_percent": count_for(0.5),
        "cells_for_80_percent": count_for(0.8),
        "radial_bands_about_contact_r_0p5": radial_bands,
        "top_10": top,
    }


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
    if coarse_report["time"] != 0.05 or fine_report["time"] != 0.05:
        raise SystemExit("expected matching t=.05 checkpoints")
    if not np.array_equal(coarse["ids"], fine["ids"]):
        raise SystemExit("cell IDs differ")
    if not np.array_equal(coarse["actual_mass"], fine["actual_mass"]):
        raise SystemExit("fixed cell masses differ")

    mass = coarse["actual_mass"].astype(np.longdouble)
    rows = {}
    for region in ("cavity", "inner_contact"):
        if not np.array_equal(coarse[region], fine[region]):
            raise SystemExit(f"region mask differs: {region}")
        mask = coarse[region]
        pressure_delta = fine["relative_pressure_error"] - coarse["relative_pressure_error"]
        density_delta = (fine["relative_density_reference_error"]
                         - coarse["relative_density_reference_error"])
        rows[region] = {
            "pressure_error_difference": concentration(
                coarse["ids"], coarse["radius"], fine["radius"],
                coarse["relative_pressure_error"], fine["relative_pressure_error"],
                pressure_delta, mass, mask),
            "density_error_difference": concentration(
                coarse["ids"], coarse["radius"], fine["radius"],
                coarse["relative_density_reference_error"],
                fine["relative_density_reference_error"], density_delta, mass, mask),
        }

    result = {
        "status": "completed_saved_t05_error_concentration_analysis",
        "time": 0.05,
        "regions": rows,
        "source_coarse": coarse_source,
        "source_fine": fine_source,
        "physical_providers": 0,
        "accepted_time_steps": 0,
        "limits": (
            "Concentration of coarse/fine differences in errors relative to each "
            "checkpoint's own moving-cell steady reference; not a causal source "
            "decomposition or convergence-order measurement."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "status": result["status"],
        "time": result["time"],
        "regions": {
            name: {
                metric: {
                    "weighted_RMS": row[metric]["weighted_RMS"],
                    "cells_for_50_percent": row[metric]["cells_for_50_percent"],
                    "cells_for_80_percent": row[metric]["cells_for_80_percent"],
                }
                for metric in ("pressure_error_difference", "density_error_difference")
            }
            for name, row in rows.items()
        },
    }, sort_keys=True))


if __name__ == "__main__":
    main()

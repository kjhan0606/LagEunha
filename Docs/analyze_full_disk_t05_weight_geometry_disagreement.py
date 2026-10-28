"""Compare saved t=.05 Laguerre weights, cell areas and physical errors.

Read-only analysis. The per-cell state slot 2 is the additive power weight W;
the common mode of W is a gauge and is removed before comparing branches.
This report measures associations only; it does not claim causality.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from compare_full_moment_checkpoint_errors import load


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while block := stream.read(1024 * 1024):
            h.update(block)
    return h.hexdigest()


def weighted_rms(values, weights, mask):
    v = np.asarray(values, dtype=np.longdouble)[mask]
    w = np.asarray(weights, dtype=np.longdouble)[mask]
    return float(np.sqrt(np.sum(w * v * v) / np.sum(w)))


def weighted_corr(a, b, weights, mask):
    x = np.asarray(a, dtype=np.longdouble)[mask]
    y = np.asarray(b, dtype=np.longdouble)[mask]
    w = np.asarray(weights, dtype=np.longdouble)[mask]
    w /= w.sum()
    x -= np.sum(w * x)
    y -= np.sum(w * y)
    xx = np.sum(w * x * x)
    yy = np.sum(w * y * y)
    if xx == 0 or yy == 0:
        return None
    return float(np.sum(w * x * y) / np.sqrt(xx * yy))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=Path(
        "/gpfs/kjhan/LagForce/Codex_full_disk_continuous_two_dt_to_t1_20261019"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"refusing to overwrite {args.output}")

    branches = {}
    specs = (("coarse", "coarse_0010"), ("fine", "fine_0020"))
    for name, step in specs:
        folder = args.run_root / name / "checkpoints" / step
        checkpoint = folder / "accepted.npz"
        with np.load(checkpoint) as data:
            ids = data["ids"].copy()
            mass = data["mass"].copy()
            state = data["state"].copy()
            seed = data["weight_rate_seed"].copy()
        n = len(ids)
        if state.size != 20 * n + 4 or seed.size != n:
            raise SystemExit(f"unexpected native-state layout: {checkpoint}")
        fields, report, source = load(folder / "physics")
        if report["time"] != 0.05 or not np.array_equal(ids, fields["ids"]):
            raise SystemExit(f"unexpected checkpoint identity: {checkpoint}")
        branches[name] = {
            "ids": ids,
            "mass": mass,
            "weight": state[:20 * n].reshape(n, 20)[:, 2],
            "weight_rate_seed": seed,
            "fields": fields,
            "source": source,
            "accepted_checkpoint": str(checkpoint.resolve()),
            "accepted_checkpoint_sha256": sha256(checkpoint),
        }

    coarse = branches["coarse"]
    fine = branches["fine"]
    if not np.array_equal(coarse["ids"], fine["ids"]):
        raise SystemExit("cell IDs differ")
    if not np.array_equal(coarse["mass"], fine["mass"]):
        raise SystemExit("fixed masses differ")
    a = coarse["fields"]["area"]
    b = fine["fields"]["area"]
    mean_area = 0.5 * (a + b)
    delta_weight_raw = fine["weight"] - coarse["weight"]
    gauge_shift = float(np.mean(delta_weight_raw))
    delta_weight = delta_weight_raw - gauge_shift
    delta_weight_area = delta_weight / mean_area
    delta_seed_raw = fine["weight_rate_seed"] - coarse["weight_rate_seed"]
    seed_gauge_shift = float(np.mean(delta_seed_raw))
    delta_seed_area = (delta_seed_raw - seed_gauge_shift) / mean_area
    delta_area = (b - a) / mean_area
    delta_actual_density = (
        (fine["fields"]["actual_density"] - coarse["fields"]["actual_density"])
        / (0.5 * (fine["fields"]["actual_density"]
                  + coarse["fields"]["actual_density"]))
    )
    delta_density_error = (
        fine["fields"]["relative_density_reference_error"]
        - coarse["fields"]["relative_density_reference_error"]
    )
    delta_pressure_error = (
        fine["fields"]["relative_pressure_error"]
        - coarse["fields"]["relative_pressure_error"]
    )
    weights = coarse["mass"].astype(np.longdouble)

    regions = {}
    for name in ("all", "cavity", "inner_contact", "disk_core"):
        mask = coarse["fields"][name]
        if not np.array_equal(mask, fine["fields"][name]):
            raise SystemExit(f"region mask differs: {name}")
        order = np.flatnonzero(mask)
        order = order[np.argsort(np.abs(delta_density_error[order]))[::-1]]
        top = []
        for i in order[:10]:
            top.append({
                "id": int(coarse["ids"][i]),
                "radius_fine": float(fine["fields"]["radius"][i]),
                "delta_weight_over_mean_area": float(delta_weight_area[i]),
                "relative_area_difference": float(delta_area[i]),
                "actual_density_relative_difference": float(delta_actual_density[i]),
                "relative_density_error_difference": float(delta_density_error[i]),
                "relative_pressure_error_difference": float(delta_pressure_error[i]),
                "gauge_removed_weight_rate_seed_difference_over_area": float(delta_seed_area[i]),
            })
        regions[name] = {
            "cells": int(mask.sum()),
            "delta_weight_over_mean_area_mass_weighted_RMS": weighted_rms(
                delta_weight_area, weights, mask),
            "relative_area_difference_mass_weighted_RMS": weighted_rms(
                delta_area, weights, mask),
            "actual_density_relative_difference_mass_weighted_RMS": weighted_rms(
                delta_actual_density, weights, mask),
            "relative_density_error_difference_mass_weighted_RMS": weighted_rms(
                delta_density_error, weights, mask),
            "relative_pressure_error_difference_mass_weighted_RMS": weighted_rms(
                delta_pressure_error, weights, mask),
            "gauge_removed_weight_rate_seed_difference_over_area_mass_weighted_RMS": weighted_rms(
                delta_seed_area, weights, mask),
            "corr_weight_difference_area_difference": weighted_corr(
                delta_weight_area, delta_area, weights, mask),
            "corr_weight_difference_relative_density_error_difference": weighted_corr(
                delta_weight_area, delta_density_error, weights, mask),
            "corr_weight_difference_relative_pressure_error_difference": weighted_corr(
                delta_weight_area, delta_pressure_error, weights, mask),
            "top_density_error_difference_cells": top,
        }

    actual_density_closure = max(
        float(np.max(np.abs(coarse["fields"]["actual_density"]
                            - coarse["mass"] / a)
                       / coarse["fields"]["actual_density"])),
        float(np.max(np.abs(fine["fields"]["actual_density"]
                            - fine["mass"] / b)
                       / fine["fields"]["actual_density"])),
    )
    result = {
        "status": "completed_saved_t05_weight_geometry_disagreement",
        "time": 0.05,
        "weight_state_slot": 2,
        "weight_unit": "length^2 additive Laguerre power weight",
        "common_weight_gauge_removed_by_unweighted_global_mean": gauge_shift,
        "common_weight_rate_seed_gauge_removed_by_unweighted_global_mean": seed_gauge_shift,
        "mass_arrays_bitwise_identical": True,
        "actual_density_equals_mass_over_area_max_relative_residual": actual_density_closure,
        "regions": regions,
        "source_coarse": {
            "checkpoint": coarse["accepted_checkpoint"],
            "checkpoint_sha256": coarse["accepted_checkpoint_sha256"],
            **coarse["source"],
        },
        "source_fine": {
            "checkpoint": fine["accepted_checkpoint"],
            "checkpoint_sha256": fine["accepted_checkpoint_sha256"],
            **fine["source"],
        },
        "physical_providers": 0,
        "accepted_time_steps": 0,
        "limits": (
            "Weight-area/error associations at one pair of coarse/fine endpoints "
            "are not a causal decomposition. Weight rates are stored seeds; the "
            "difference does not prove they caused the density/pressure error."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "status": result["status"],
        "time": result["time"],
        "gauge_shift_removed": gauge_shift,
        "actual_density_mass_area_max_residual": actual_density_closure,
        "regions": {
            name: {
                key: row[key]
                for key in (
                    "delta_weight_over_mean_area_mass_weighted_RMS",
                    "relative_area_difference_mass_weighted_RMS",
                    "relative_density_error_difference_mass_weighted_RMS",
                    "relative_pressure_error_difference_mass_weighted_RMS",
                    "corr_weight_difference_area_difference",
                    "corr_weight_difference_relative_density_error_difference",
                    "corr_weight_difference_relative_pressure_error_difference",
                )
            }
            for name, row in regions.items()
        },
    }, sort_keys=True))


if __name__ == "__main__":
    main()

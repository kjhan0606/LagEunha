"""Condition saved spin/thermodynamic errors on broad spatial trends.

This is a read-only follow-up to analyze_saved_rotation_error_cospatiality.py.
It reports mass-weighted RMS errors and partial correlations after removing a
low-order radial trend plus first/second azimuthal harmonics. These diagnostics
test whether the observed association is only a smooth spatial pattern; they do
not establish causation or long-time accuracy.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from compare_full_moment_checkpoint_errors import load


CASES = (
    ("t0025_coarse", "coarse", "coarse_0005"),
    ("t0025_fine", "fine", "fine_0010"),
    ("t005_coarse", "coarse", "coarse_0010"),
    ("t005_fine", "fine", "fine_0020"),
)
REGIONS = ("inner_contact", "cavity", "disk_core")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def weighted_corr(a, b, weight):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    w = np.asarray(weight, dtype=np.float64)
    w = w / np.sum(w)
    a = a - np.sum(w * a)
    b = b - np.sum(w * b)
    aa = np.sum(w * a * a)
    bb = np.sum(w * b * b)
    return None if aa == 0 or bb == 0 else float(np.sum(w*a*b)/np.sqrt(aa*bb))


def detrend(values, design, weight):
    values = np.asarray(values, dtype=np.float64)
    design = np.asarray(design, dtype=np.float64)
    weight = np.asarray(weight, dtype=np.float64)
    sw = np.sqrt(weight / np.sum(weight))
    coeff, *_ = np.linalg.lstsq(design * sw[:, None], values * sw, rcond=None)
    return values - design @ coeff


def weighted_rms(values, weight):
    values = np.asarray(values, dtype=np.float64)
    weight = np.asarray(weight, dtype=np.float64)
    return float(np.sqrt(np.average(values*values, weights=weight)))


def summarize_region(ids, mass, mask, centroid, spin, density, thermal, pressure):
    idx = np.flatnonzero(mask)
    if not len(idx):
        return {"cells": 0}
    w = np.asarray(mass[idx], dtype=np.float64)
    xy = np.asarray(centroid[idx], dtype=np.float64)
    radius = np.linalg.norm(xy, axis=1)
    theta = np.arctan2(xy[:, 1], xy[:, 0])
    rbar = np.average(radius, weights=w)
    rsigma = np.sqrt(np.average((radius-rbar)**2, weights=w))
    rr = (radius-rbar) / max(rsigma, np.finfo(float).tiny)
    design = np.column_stack((
        np.ones(len(idx)), rr, rr*rr,
        np.sin(theta), np.cos(theta), np.sin(2*theta), np.cos(2*theta),
    ))
    local_spin = np.asarray(spin[idx], dtype=np.float64)
    errors = {
        "relative_density": np.asarray(density[idx], dtype=np.float64),
        "relative_specific_thermal": np.asarray(thermal[idx], dtype=np.float64),
        "relative_pressure": np.asarray(pressure[idx], dtype=np.float64),
    }
    correlations = {}
    for name, error in errors.items():
        correlations[name] = {
            "raw": weighted_corr(local_spin, error, w),
            "after_radial_and_azimuthal_detrending": weighted_corr(
                detrend(local_spin, design, w), detrend(error, design, w), w),
        }
    return {
        "cells": int(len(idx)),
        "mass_weighted_rms": {
            "intrinsic_spin_error": weighted_rms(local_spin, w),
            **{name: weighted_rms(error, w) for name, error in errors.items()},
        },
        "density_thermal_cancellation": {
            "relative_density_vs_specific_thermal_error_correlation": weighted_corr(
                errors["relative_density"], errors["relative_specific_thermal"], w),
            "linear_sum_rms": weighted_rms(
                errors["relative_density"] + errors["relative_specific_thermal"], w),
        },
        "spin_error_correlations": correlations,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=Path(
        "/gpfs/kjhan/LagForce/Codex_full_disk_continuous_two_dt_to_t1_20261019"))
    parser.add_argument("--rotation-root", type=Path, default=Path("Docs"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"refusing to overwrite {args.output}")

    cases = {}
    for tag, branch, checkpoint in CASES:
        physics = args.run_root / branch / "checkpoints" / checkpoint / "physics"
        fields, report, physical_source = load(physics)
        rotation_dir = args.rotation_root / f"canonical_saved_cell_rotation_{tag}_2026-10-19"
        rotation_path = rotation_dir / "rotation.npz"
        with np.load(rotation_path) as z:
            ids = z["ids"].copy()
            spin = np.asarray(z["actual_ell"] - z["reference_ell"], dtype=np.float64)
        if not np.array_equal(ids, fields["ids"]):
            raise SystemExit(f"cell ID mismatch in {tag}")
        if np.max(np.abs(spin - (fields["actual_intrinsic_ell"]
                                 - fields["reference_intrinsic_ell"]))) > 1e-18:
            raise SystemExit(f"spin-error cross-check failed in {tag}")

        thermal_error = (fields["actual_specific_thermal"]
                         / fields["reference_specific_thermal"] - 1.0)
        pressure_identity = ((1.0 + fields["relative_density_reference_error"])
                             * (1.0 + thermal_error) - 1.0)
        pressure_identity_max_abs = float(np.max(np.abs(
            pressure_identity - fields["relative_pressure_error"])))
        if pressure_identity_max_abs > 1e-10:
            raise SystemExit(f"pressure identity cross-check failed in {tag}: "
                             f"{pressure_identity_max_abs:.3e}")

        cases[tag] = {
            "time": float(report["time"]),
            "source": {
                "physical_comparison": physical_source,
                "rotation_payload_sha256": sha256(rotation_path),
            },
            "pressure_identity_max_absolute": pressure_identity_max_abs,
            "regions": {
                name: summarize_region(
                    ids, fields["actual_mass"], fields[name], fields["actual_centroid"],
                    spin, fields["relative_density_reference_error"], thermal_error,
                    fields["relative_pressure_error"],
                ) for name in REGIONS
            },
        }

    result = {
        "status": "completed_saved_rotation_error_spatially_conditioned",
        "scope": "Four saved endpoints; descriptive cellwise analysis only.",
        "conditioning": (
            "Mass-weighted least-squares removal of constant, standardized radius, "
            "radius squared, and azimuthal Fourier modes m=1,2 on each region."
        ),
        "error_definitions": {
            "spin": "actual intrinsic ell minus steady Kepler reference ell",
            "density": "relative density error against the steady reference",
            "specific_thermal": "actual/reference specific thermal energy minus one",
            "pressure": "relative pressure error against the steady reference",
            "pressure_identity": "e_P=(1+e_rho)(1+e_u)-1 for the ideal-gas state",
        },
        "cases": cases,
        "limits": [
            "Partial correlations do not establish causal coupling.",
            "The spatial detrend removes only broad low-order radial and azimuthal structure.",
            "The four endpoints stop at t=0.05 and do not test long-term Kepler stability.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "status": result["status"],
        "cases": {
            tag: {
                region: {
                    "cells": row["cells"],
                    "rms": row["mass_weighted_rms"],
                    "density_thermal_cancellation": row["density_thermal_cancellation"],
                    "partial_corr": {
                        key: value["after_radial_and_azimuthal_detrending"]
                        for key, value in row["spin_error_correlations"].items()
                    },
                } for region, row in case["regions"].items()
            } for tag, case in cases.items()
        },
    }, sort_keys=True))


if __name__ == "__main__":
    main()

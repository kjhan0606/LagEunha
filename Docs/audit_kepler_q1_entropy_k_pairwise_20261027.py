"""Compare native q1 energy- and fixed-K RK4 runs from the same checkpoint.

The useful entropy increment begins at the shared t=.825 accepted state, not
at the t=0 IC; both arms inherit its pre-existing entropy error.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_polar_long_gate_20261027 import ordered, total_energy


GAMMA_MINUS_ONE = float(np.float32(1.666667)) - 1.0
N = 384 * 384
STRIDE = 576
K_OFFSET = 552
REGIONS = (("core", 0., .5), ("disk", .5, 2.),
           ("inner_disk", .5, .75), ("outer_disk", 1., 2.),
           ("transition", 2., 2.35), ("outer_background", 2.75, 5.))


def weighted_mean(values, weights):
    return float(np.sum(values * weights, dtype=np.longdouble)
                 / np.sum(weights, dtype=np.longdouble))


def ordered_k(path):
    """Read the stored entropy variable by cell ID, independently of EOS fields."""
    layout = np.dtype({"names": ("id_word", "K"),
                       "formats": ("<u8", "<f8"),
                       "offsets": (8, K_OFFSET), "itemsize": STRIDE})
    data = np.memmap(path, dtype=layout, mode="r", offset=16, shape=(N,))
    ids = (data["id_word"] & np.uint64((1 << 56) - 1)).astype(np.int64)
    order = np.argsort(ids)
    if not np.array_equal(ids[order], np.arange(N)):
        raise ValueError(f"Invalid cell IDs in K payload: {path}")
    return np.asarray(data["K"][order], dtype=float)


def region_audit(start, end, reference_radius, mass, lo, hi):
    mask = (reference_radius >= lo) & (reference_radius < hi)
    if not np.any(mask):
        raise ValueError(f"Empty reference region {lo}..{hi}")
    u0 = start["ie"].astype(float)
    u = end["ie"].astype(float)
    v0 = start["volume"].astype(float)
    v = end["volume"].astype(float)
    entropy_increment = np.log(u / u0) + GAMMA_MINUS_ONE * np.log(v / v0)
    adiabatic_u = u0 * (v0 / v) ** GAMMA_MINUS_ONE
    excess = u - adiabatic_u
    x = end["x"].astype(float) - 3.
    y = end["y"].astype(float) - 3.
    x0 = start["x"].astype(float) - 3.
    y0 = start["y"].astype(float) - 3.
    L0 = np.sum(mass[mask] * (x0[mask] * start["vy"][mask]
                              - y0[mask] * start["vx"][mask]), dtype=np.longdouble)
    Lt = np.sum(mass[mask] * (x[mask] * end["vy"][mask]
                              - y[mask] * end["vx"][mask]), dtype=np.longdouble)
    u0_sum = float(np.sum(u0[mask], dtype=np.longdouble))
    return {
        "cells": int(np.sum(mask)),
        "mass": float(np.sum(mass[mask], dtype=np.longdouble)),
        "entropy_increment_signed_mass_weighted": weighted_mean(entropy_increment[mask], mass[mask]),
        "entropy_increment_abs_mass_weighted": weighted_mean(abs(entropy_increment[mask]), mass[mask]),
        "entropy_positive_mass_fraction": float(
            np.sum(mass[mask & (entropy_increment > 0)], dtype=np.longdouble)
            / np.sum(mass[mask], dtype=np.longdouble)),
        "thermal_initial": u0_sum,
        "thermal_final": float(np.sum(u[mask], dtype=np.longdouble)),
        "adiabatic_thermal_final": float(np.sum(adiabatic_u[mask], dtype=np.longdouble)),
        "thermal_excess_signed_over_initial": float(np.sum(excess[mask], dtype=np.longdouble) / u0_sum),
        "thermal_excess_abs_over_initial": float(np.sum(abs(excess[mask]), dtype=np.longdouble) / u0_sum),
        "angular_momentum_relative": float(Lt / L0 - 1.) if L0 != 0 else None,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((args.root / "manifest.json").read_text())
    reference, t0 = ordered(Path(manifest["initial"]))
    start, ts = ordered(Path(manifest["parent_checkpoint"]))
    if t0 != 0 or abs(ts - .825) > 1e-10:
        raise ValueError("Incorrect shared initial and start checkpoint")
    mass = reference["mass"].astype(float)
    if not np.array_equal(reference["mass"], start["mass"]):
        raise ValueError("Shared start changed initial cell mass")
    radius = np.hypot(reference["x"].astype(float) - 3.,
                      reference["y"].astype(float) - 3.)
    disk_mask = (radius >= .5) & (radius < 2.)
    cases = {}
    endpoints = {}
    for label, mode in (("energy", 0), ("entropy", 1)):
        run = args.root / label
        if manifest["runs"][label]["entropy_mode"] != mode:
            raise ValueError("Staged mode mismatch")
        receipt = json.loads((run / "audit_t0845.json").read_text())
        radial_path = run / "audit_radial_t0845.json"
        radial = json.loads(radial_path.read_text()) if radial_path.exists() else None
        if abs(receipt["time"] - .845) > 1e-10:
            raise ValueError(f"Native endpoint time mismatch in {label}")
        endpoint = Path(receipt["endpoint"])
        if endpoint.parent.resolve() != run.resolve():
            raise ValueError("Receipt endpoint is outside staged run")
        end, te = ordered(endpoint)
        if abs(te - .845) > 1e-10 or not np.array_equal(start["mass"], end["mass"]):
            raise ValueError("Incompatible endpoint time or mass")
        if np.any(end["w2"] != 0) or np.any(end["ie"] <= 0) or np.any(end["volume"] <= 0):
            raise ValueError("W0 or positive state contract failed")
        log = (run / "log").read_text(errors="replace")
        if "EXIT:0" not in log:
            raise ValueError("Native executable did not exit successfully")
        banner = "[ENTROPY_K_INIT] rank=0 initialized=147456 retained=0"
        if (banner in log) != (mode == 1):
            raise ValueError("Entropy-K initialization banner mismatch")
        endpoints[label] = end
        stored_k = ordered_k(endpoint)
        if not np.isfinite(stored_k).all():
            raise ValueError(f"Nonfinite stored K in {label}")
        endpoint_u = end["ie"].astype(float)
        endpoint_v = end["volume"].astype(float)
        endpoint_p = end["pressure"].astype(float)
        endpoint_rho = end["mass"].astype(float) / endpoint_v
        thermal_eos_error = float(np.max(abs(
            endpoint_u / (endpoint_p * endpoint_v / GAMMA_MINUS_ONE) - 1.)))
        if not np.isfinite(thermal_eos_error) or thermal_eos_error > 1e-5:
            raise ValueError(f"U/P EOS mismatch in {label}: {thermal_eos_error}")
        k_info = {"stored_nonzero": int(np.count_nonzero(stored_k))}
        if mode == 1:
            u0 = start["ie"].astype(float)
            v0 = start["volume"].astype(float)
            rho0 = mass / v0
            k_expected = GAMMA_MINUS_ONE * (u0 / v0) / rho0 ** (GAMMA_MINUS_ONE + 1.)
            relative_k = stored_k / k_expected - 1.
            if np.any(stored_k <= 0):
                raise ValueError("Entropy-mode K is not positive in every cell")
            k_info.update({
                "relative_K_from_accepted_state_abs_max": float(np.max(abs(relative_k))),
                "relative_K_from_accepted_state_disk_abs_mass_weighted":
                    weighted_mean(abs(relative_k[disk_mask]), mass[disk_mask]),
            })
            pressure_k_error = float(np.max(abs(
                endpoint_p / (stored_k * endpoint_rho ** (GAMMA_MINUS_ONE + 1.)) - 1.)))
            k_info["pressure_from_stored_K_relative_abs_max"] = pressure_k_error
            if not np.isfinite(pressure_k_error) or pressure_k_error > 1e-5:
                raise ValueError(f"P/K EOS mismatch in {label}: {pressure_k_error}")
        elif np.count_nonzero(stored_k):
            raise ValueError("Energy-mode control unexpectedly changed stored K")
        cases[label] = {
            "endpoint": str(endpoint), "time": te,
            "mode": mode,
            "total_energy_relative_from_start": (
                total_energy(end) - total_energy(start)) / abs(total_energy(start)),
            "regions": {name: region_audit(start, end, radius, mass, lo, hi)
                        for name, lo, hi in REGIONS},
            "native_gate_checks": receipt["checks"],
            "native_gate_pass": bool(all(receipt["checks"].values())),
            "radial_gate_checks": radial["checks"] if radial else None,
            "radial_gate_pass": bool(all(radial["checks"].values())) if radial else None,
            "final_ledger": receipt["rk4_final_ledger"],
            "stored_K": k_info,
            "thermal_U_from_pressure_EOS_relative_abs_max": thermal_eos_error,
        }
    disk = disk_mask
    old = endpoints["energy"]
    new = endpoints["entropy"]
    result = {
        "manifest": str(args.root / "manifest.json"),
        "start_time": ts, "end_time": .845,
        "cases": cases,
        "same_time_disk_difference": {
            "position_abs_mass_weighted": weighted_mean(
                np.hypot(new["x"][disk] - old["x"][disk],
                         new["y"][disk] - old["y"][disk]), mass[disk]),
            "velocity_abs_mass_weighted": weighted_mean(
                np.hypot(new["vx"][disk] - old["vx"][disk],
                         new["vy"][disk] - old["vy"][disk]), mass[disk]),
            "density_relative_abs_mass_weighted": weighted_mean(
                abs(new["den"][disk] / old["den"][disk] - 1.), mass[disk]),
        },
        "scope": "Matched 0.02-time native RK4 attribution comparison from the same accepted q1 checkpoint; failed physical gates remain failures and neither arm establishes t>6 stability or a speedup.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({
        "energy_disk": cases["energy"]["regions"]["disk"],
        "entropy_disk": cases["entropy"]["regions"]["disk"],
        "energy_drift": cases["energy"]["total_energy_relative_from_start"],
        "entropy_drift": cases["entropy"]["total_energy_relative_from_start"],
        "same_time_disk_difference": result["same_time_disk_difference"],
    }, indent=2))


if __name__ == "__main__":
    main()

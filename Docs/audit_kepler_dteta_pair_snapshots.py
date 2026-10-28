#!/usr/bin/env python3
"""Compare capped/no-cap Kepler snapshots to one common W0 initial state."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

SNAPSHOT_TOOLS = Path(
    "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020"
)
sys.path.insert(0, str(SNAPSHOT_TOOLS))
from audit_snapshot import read_snapshot  # noqa: E402


def load(path):
    particles, meta = read_snapshot(path)
    order = np.argsort(particles["u4if"])
    return particles[order], meta


def state(path, center, eps, bins):
    particles, meta = load(path)
    for field in ("x", "y", "vx", "vy", "mass", "volume", "ie",
                  "den", "pressure", "w2"):
        if not np.isfinite(np.asarray(particles[field], dtype=np.float64)).all():
            raise ValueError(f"{path}: nonfinite values in {field}")
    mass = np.asarray(particles["mass"], dtype=np.float64)
    if np.any(mass <= 0):
        raise ValueError(f"{path}: nonpositive cell mass")
    x = np.asarray(particles["x"], dtype=np.float64) - center[0]
    y = np.asarray(particles["y"], dtype=np.float64) - center[1]
    vx = np.asarray(particles["vx"], dtype=np.float64)
    vy = np.asarray(particles["vy"], dtype=np.float64)
    radius = np.hypot(x, y)
    safe_radius = np.maximum(radius, 1e-30)
    radial_velocity = (x * vx + y * vy) / safe_radius
    azimuthal_velocity = (x * vy - y * vx) / safe_radius
    circular_velocity = radius * (radius * radius + eps * eps) ** (-0.75)
    angular_momentum = x * vy - y * vx
    volume = np.asarray(particles["volume"], dtype=np.float64)
    internal_energy = np.asarray(particles["ie"], dtype=np.float64)
    total_energy = float(
        internal_energy.sum()
        + 0.5 * np.sum(mass * (vx * vx + vy * vy))
        - np.sum(mass / np.sqrt(radius * radius + eps * eps))
    )

    density_profile = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        selected = (radius >= lo) & (radius < hi)
        if np.any(selected) and np.any(volume[selected] <= 0):
            raise ValueError(f"{path}: nonpositive volume in density bin {lo}:{hi}")
        density_profile.append(
            float(mass[selected].sum() / volume[selected].sum())
            if np.any(selected)
            else None
        )
    values = {
        "ids": np.asarray(particles["u4if"]),
        "mass": mass,
        "w2": np.asarray(particles["w2"], dtype=np.float64),
        "ie": internal_energy,
        "volume": volume,
        "radius": radius,
        "radial_velocity": radial_velocity,
        "azimuthal_velocity": azimuthal_velocity,
        "circular_velocity": circular_velocity,
        "angular_momentum": angular_momentum,
        "density_profile": density_profile,
        "total_energy": total_energy,
    }
    return values, meta


def metrics(label, path, initial, initial_radius, active, center, eps, bins):
    values, meta = state(path, center, eps, bins)
    if not np.array_equal(values["ids"], initial["ids"]):
        raise ValueError(f"{label}: particle IDs differ from initial snapshot")
    if not np.array_equal(values["mass"], initial["mass"]):
        raise ValueError(f"{label}: per-ID fixed-mass invariant failed")
    if np.any(values["w2"] != 0):
        raise ValueError(f"{label}: endpoint is not W=0")
    if np.any(values["ie"] <= 0) or np.any(values["volume"] <= 0):
        raise ValueError(f"{label}: nonpositive internal energy or volume")

    profile = np.asarray(values["density_profile"], dtype=object)
    initial_profile = np.asarray(initial["density_profile"], dtype=object)
    valid = np.array(
        [a is not None and b is not None for a, b in zip(profile, initial_profile)]
    )
    rho = profile.astype(float)[valid]
    rho0 = initial_profile.astype(float)[valid]
    lz0_disk = np.sum(initial["mass"][active] * initial["angular_momentum"][active])
    lzt_disk = np.sum(values["mass"][active] * values["angular_momentum"][active])
    lz0_all = np.sum(initial["mass"] * initial["angular_momentum"])
    lzt_all = np.sum(values["mass"] * values["angular_momentum"])
    return {
        "snapshot": str(path),
        "time": float(meta["t"]),
        "cell_count": int(len(values["ids"])),
        "mass_exact_by_id": True,
        "W0": True,
        "active_count": int(active.sum()),
        "active_radius_error_rms": float(
            np.sqrt(np.mean(((values["radius"][active] - initial_radius[active])
                             / initial_radius[active]) ** 2))
        ),
        "active_vr_over_vc_rms": float(
            np.sqrt(np.mean((values["radial_velocity"][active]
                             / values["circular_velocity"][active]) ** 2))
        ),
        "active_abs_vphi_error_over_vc_mean": float(
            np.mean(np.abs(values["azimuthal_velocity"][active]
                           - values["circular_velocity"][active])
                    / values["circular_velocity"][active])
        ),
        "active_Lz_relative_change": float(lzt_disk / lz0_disk - 1.0),
        "global_Lz_relative_change": float(lzt_all / lz0_all - 1.0),
        "relative_total_energy_change_from_initial": float(
            (values["total_energy"] - initial["total_energy"])
            / abs(initial["total_energy"])
        ),
        "annular_density_profile_mean_abs_relative_change": float(
            np.mean(np.abs((rho - rho0) / rho0))
        ) if len(rho) else None,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--initial", type=Path, required=True)
    parser.add_argument("--capped", type=Path, required=True)
    parser.add_argument("--uncapped", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--center", type=float, nargs=2, default=(3.0, 3.0))
    parser.add_argument("--eps", type=float, default=0.01175)
    parser.add_argument("--r-inner", type=float, default=0.5)
    parser.add_argument("--r-outer", type=float, default=2.0)
    parser.add_argument("--bins", type=int, default=15)
    args = parser.parse_args()

    edges = np.linspace(args.r_inner, args.r_outer, args.bins + 1)
    initial, initial_meta = state(args.initial, args.center, args.eps, edges)
    initial_radius = initial["radius"]
    active = (initial_radius >= args.r_inner) & (initial_radius <= args.r_outer)
    result = {
        "method": "common-initial fixed-mass W0 Kepler kinematic comparison",
        "initial_snapshot": str(args.initial),
        "initial_time": float(initial_meta["t"]),
        "center": list(args.center),
        "softening": args.eps,
        "active_annulus": [args.r_inner, args.r_outer],
        "density_bins": edges.tolist(),
        "capped": metrics("capped", args.capped, initial, initial_radius,
                           active, args.center, args.eps, edges),
        "uncapped": metrics("uncapped", args.uncapped, initial, initial_radius,
                             active, args.center, args.eps, edges),
    }
    result["endpoint_time_difference"] = (
        result["uncapped"]["time"] - result["capped"]["time"]
    )
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

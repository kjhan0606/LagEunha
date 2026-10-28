"""Compare complete equal-time W=0 Kepler checkpoints with dteta .06/.03."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_q2_weighted_checkpoint_pair_20261027 import measure
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N = 55351
GAMMA_MINUS_ONE = float(np.float32(1.666667)) - 1.0


def mass_mean_abs(values, mass):
    return float(np.sum(np.abs(values) * mass, dtype=np.longdouble) /
                 np.sum(mass, dtype=np.longdouble))


def orbit_residual(p0, p, t):
    x0 = p0["x"].astype(float) - 3.0
    y0 = p0["y"].astype(float) - 3.0
    r0 = np.hypot(x0, y0)
    phase = t * (r0 * r0 + 0.01175 ** 2) ** (-0.75)
    exact_x = x0 * np.cos(phase) - y0 * np.sin(phase)
    exact_y = x0 * np.sin(phase) + y0 * np.cos(phase)
    return np.column_stack((p["x"].astype(float) - 3.0 - exact_x,
                            p["y"].astype(float) - 3.0 - exact_y))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--standard", type=Path, required=True)
    ap.add_argument("--half", type=Path, required=True)
    ap.add_argument("--initial", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    runs = {"dteta_0.06": args.standard.resolve(), "dteta_0.03": args.half.resolve()}
    p0, _ = ordered(args.initial, N)
    initial_sha = digest(args.initial)
    if abs(snapshot_time(args.initial, N)) > 1e-12:
        raise ValueError("Reference must be the t=0 IC")
    start = {}
    state = {}
    time = {}
    for label, path in runs.items():
        restart = path.parent / "khout.000001.dat"
        binary = path.parent / "eunha2.exe"
        pr, _ = ordered(restart, N)
        p, _ = ordered(path, N)
        if not (np.array_equal(pr["mass"], p0["mass"]) and
                np.array_equal(p["mass"], p0["mass"]) and
                np.all(p["w2"] == 0) and
                np.all(p["volume"] > 0) and
                np.all(p["den"] > 0) and
                np.all(p["ie"] > 0)):
            raise ValueError(f"Invalid W=0 fixed-mass checkpoint: {label}")
        if abs(snapshot_time(restart, N) - 6.5) > 1e-10:
            raise ValueError(f"Unexpected restart time: {label}")
        start[label] = (digest(restart), digest(binary))
        state[label] = p
        time[label] = snapshot_time(path, N)

    first, second = tuple(runs)
    if start[first] != start[second]:
        raise ValueError("The t=6.5 restart or native executable differs between timesteps")
    if abs(time[first] - time[second]) > 1e-6:
        raise ValueError("Snapshots are not at the same physical time")

    radius0 = np.hypot(p0["x"].astype(float) - 3, p0["y"].astype(float) - 3)
    bands = {}
    for name, lo, hi in (("inner", 0.5, 0.75), ("disk", 0.5, 2.0)):
        select = (radius0 >= lo) & (radius0 < hi)
        a = measure(p0, state[first], time[first], select)
        b = measure(p0, state[second], time[second], select)
        mass = p0["mass"][select].astype(float)
        den0 = p0["den"][select].astype(float)
        pa, pb = state[first], state[second]
        density_delta = (pb["den"][select].astype(float) -
                         pa["den"][select].astype(float)) / den0
        velocity_delta = np.hypot(pb["vx"][select].astype(float) -
                                  pa["vx"][select].astype(float),
                                  pb["vy"][select].astype(float) -
                                  pa["vy"][select].astype(float))
        entropy_delta = np.log(pb["ie"][select].astype(float) /
                               pa["ie"][select].astype(float)) + (
            GAMMA_MINUS_ONE *
            np.log(pb["volume"][select].astype(float) /
                   pa["volume"][select].astype(float)))
        position_delta = np.linalg.norm(
            orbit_residual(p0[select], pb[select], time[second]) -
            orbit_residual(p0[select], pa[select], time[first]), axis=1)
        bands[name] = {
            "radius_range": [lo, hi], first: a, second: b,
            "half_over_standard": {
                key: b[key] / a[key]
                for key in ("density_abs_mass_mean", "entropy_log_abs_mass_mean",
                            "orbit_position_error_abs_mass_mean")
            },
            "between_timesteps": {
                "density_relative_abs_mass_mean": mass_mean_abs(density_delta, mass),
                "density_relative_abs_p99": float(np.quantile(np.abs(density_delta), 0.99)),
                "velocity_vector_abs_mass_mean": mass_mean_abs(velocity_delta, mass),
                "entropy_log_abs_mass_mean": mass_mean_abs(entropy_delta, mass),
                "orbit_residual_vector_abs_mass_mean": mass_mean_abs(position_delta, mass),
            },
        }

    result = {
        "runs": {
            label: {"checkpoint": str(path), "checkpoint_sha256": digest(path),
                    "time": time[label]}
            for label, path in runs.items()
        },
        "initial_sha256": initial_sha,
        "common_t65_restart_sha256": start[first][0],
        "common_executable_sha256": start[first][1],
        "time_difference": abs(time[first] - time[second]),
        "bands": bands,
        "scope": "Read-only native completed-checkpoint comparison from the accepted t=6.5 restart. No endpoint audit or convergence claim from one dt pair.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

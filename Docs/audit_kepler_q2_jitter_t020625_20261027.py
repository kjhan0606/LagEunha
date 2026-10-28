"""Compare native q2 and 5% site-jitter runs at one identical saved time."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_polar_checkpoint_20261026 import (
    CX, CY, GAMMA_MINUS_ONE, analyze, load_ordered, weighted_mean,
)
from prepare_kepler_q2_jitter_to020625_20261027 import TARGET
from prepare_kepler_polar_long_restart_20261027 import (
    digest, require_zero_injection, snapshot_time,
)


def disk_metrics(initial, endpoint):
    p0, _ = load_ordered(initial)
    p, _ = load_ordered(endpoint)
    r0 = np.hypot(p0["x"] - CX, p0["y"] - CY)
    disk = (r0 >= .5) & (r0 < 2.)
    mass = p0["mass"][disk].astype(float)
    rho_change = p["den"][disk].astype(float) / p0["den"][disk].astype(float) - 1
    entropy_change = (p["ie"][disk].astype(float) /
                      (p0["ie"][disk].astype(float) *
                       (p0["volume"][disk].astype(float) /
                        p["volume"][disk].astype(float)) ** GAMMA_MINUS_ONE) - 1)
    x0, y0 = p0["x"][disk] - CX, p0["y"][disk] - CY
    x1, y1 = p["x"][disk] - CX, p["y"][disk] - CY
    vr = (x1 * p["vx"][disk] + y1 * p["vy"][disk]) / np.hypot(x1, y1)
    l0 = np.sum(mass * (x0 * p0["vy"][disk] - y0 * p0["vx"][disk]), dtype=np.longdouble)
    l1 = np.sum(mass * (x1 * p["vy"][disk] - y1 * p["vx"][disk]), dtype=np.longdouble)
    return {"cells": int(disk.sum()),
            "mass": float(np.sum(mass, dtype=np.longdouble)),
            "density_change_signed_mass_weighted": weighted_mean(rho_change, mass),
            "density_change_abs_mass_weighted": weighted_mean(abs(rho_change), mass),
            "entropy_change_signed_mass_weighted": weighted_mean(entropy_change, mass),
            "entropy_change_abs_mass_weighted": weighted_mean(abs(entropy_change), mass),
            "radial_velocity_abs_mass_weighted": weighted_mean(abs(vr), mass),
            "angular_momentum_relative": float(l1 / l0 - 1),
            "floor_free_positive_state": bool(np.all(p["ie"] > 0) and np.all(p["volume"] > 0)),
            "mass_exact": bool(np.array_equal(p0["mass"], p["mass"])),
            "W0_exact": bool(np.all(p0["w2"] == 0) and np.all(p["w2"] == 0))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("jitter_preflight", type=Path)
    parser.add_argument("jitter_run", type=Path)
    parser.add_argument("q2_initial", type=Path)
    parser.add_argument("q2_endpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    jitter_initial = args.jitter_preflight / "khout.000000.dat"
    snapshots = sorted(args.jitter_run.glob("khout.*.dat"))
    if len(snapshots) != 2:
        raise ValueError("Expected only one restart checkpoint and one endpoint")
    jitter_endpoint = snapshots[-1]
    times = (snapshot_time(jitter_initial), snapshot_time(jitter_endpoint),
             snapshot_time(args.q2_initial), snapshot_time(args.q2_endpoint))
    if (abs(times[0]) > 1e-12 or abs(times[2]) > 1e-12 or
            abs(times[1] - TARGET) > 1e-10 or abs(times[3] - TARGET) > 1e-10):
        raise ValueError(f"Unmatched native snapshot times: {times}")
    require_zero_injection(args.jitter_run / "log")
    cases = {}
    for name, initial, endpoint in (
        ("q2", args.q2_initial, args.q2_endpoint),
        ("q2_site_jitter_005", jitter_initial, jitter_endpoint),
    ):
        detail = analyze(initial, endpoint)
        if not all(detail["checks"].values()):
            raise ValueError(f"{name} failed physical state gate")
        disk = disk_metrics(initial, endpoint)
        if not all(disk[k] for k in ("floor_free_positive_state", "mass_exact", "W0_exact")):
            raise ValueError(f"{name} failed fixed-mass W0 positive state gate")
        cases[name] = {"initial": str(initial), "endpoint": str(endpoint),
                       "endpoint_sha256": digest(endpoint),
                       "disk": disk,
                       "early_geometry_and_radial_bins": detail["bins"]}
    base = cases["q2"]["disk"]
    jitter = cases["q2_site_jitter_005"]["disk"]
    ratios = {key: jitter[key] / base[key] if base[key] != 0 else None
              for key in ("density_change_abs_mass_weighted",
                          "entropy_change_abs_mass_weighted",
                          "radial_velocity_abs_mass_weighted")}
    report = {"matched_time": TARGET, "cases": cases,
              "jitter_over_q2_abs_metric_ratios": ratios,
              "scope": "Matched native W0 RK4+CUDA early-time IC sensitivity control. This does not establish t>6 or t=600 stability."}
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"matched_time": TARGET, "disk": {name: case["disk"] for name, case in cases.items()},
                      "ratios": ratios}, indent=2))


if __name__ == "__main__":
    main()

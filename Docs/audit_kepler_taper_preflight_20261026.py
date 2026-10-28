"""Verify that a native sparse Kepler preflight used the outer taper IC."""

import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot


def main():
    run = Path.cwd()
    initial = run / "khout.000000.dat"
    p, meta = read_snapshot(initial)
    if float(meta["t"]) != 0.0 or (float(meta["Lx"]), float(meta["Ly"])) != (6.0, 6.0):
        raise ValueError("Expected t=0 six-by-six native Kepler snapshot")
    radius = np.hypot(p["x"].astype(float)-3.0, p["y"].astype(float)-3.0)
    inner, outer = 2.35, 2.75
    eps = 0.01175
    z = -1.0 / np.sqrt(radius**2+eps**2)
    z0 = -1.0 / np.sqrt(inner**2+eps**2)
    z1 = -1.0 / np.sqrt(outer**2+eps**2)
    u = np.clip((z-z0)/(z1-z0), 0.0, 1.0)
    speed_factor = 1.0-u**3*(10.0+u*(-15.0+6.0*u))
    speed_factor = np.where(radius<=inner, 1.0,
                            np.where(radius>=outer, 0.0, speed_factor))
    v_expected = radius * (radius**2+eps**2)**(-0.75) * speed_factor
    v_current = np.hypot(p["vx"].astype(float), p["vy"].astype(float))
    outer_mask = radius>=outer
    disk_mask = (radius>=.5)&(radius<2.0)
    if not (np.any(outer_mask) and np.any(disk_mask)):
        raise ValueError("Missing expected disk or outer-background cells")
    result = {
        "initial": str(initial), "cells": len(p),
        "outer_static_cells": int(np.sum(outer_mask)),
        "disk_cells": int(np.sum(disk_mask)),
        "outer_max_speed": float(np.max(v_current[outer_mask])),
        "speed_error_max": float(np.max(np.abs(v_current-v_expected))),
        "disk_speed_error_max": float(np.max(np.abs(v_current[disk_mask]-v_expected[disk_mask]))),
        "disk_pressure_max_relative_error": float(np.max(np.abs(p["pressure"][disk_mask]/1e-6-1.0))),
        "checks": {
            "outer_background_static": bool(np.max(v_current[outer_mask]) < 1e-10),
            "all_speeds_match_balanced_taper": bool(np.max(np.abs(v_current-v_expected)) < 2e-6),
            "disk_pressure_initial": bool(np.max(np.abs(p["pressure"][disk_mask]/1e-6-1.0)) < 1e-6),
            "W0": bool(np.all(p["w2"]==0)),
        },
        "scope": "Native t=0 runtime IC check. Continuation is needed to attribute outer perturbation or long-time disk behavior.",
    }
    with (run / "audit_taper_ic.json").open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))
    if not all(result["checks"].values()):
        raise SystemExit("Taper IC preflight gate failed")


if __name__ == "__main__":
    main()

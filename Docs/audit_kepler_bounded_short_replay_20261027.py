"""Cheap native fixed-polygon geometry replay gate from an accepted state."""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_state_consistency_20261027 import state_metrics
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


COUNT = 55351
DOMAIN_AREA = .5 * 1024 * 2.75**2 * np.sin(2 * np.pi / 1024)


def energy(p):
    x = p["x"].astype(float)-3.
    y = p["y"].astype(float)-3.
    m = p["mass"].astype(float)
    vx, vy = p["vx"].astype(float), p["vy"].astype(float)
    return float(np.sum(.5*m*(vx*vx+vy*vy)+p["ie"].astype(float)
                        -m/np.sqrt(x*x+y*y+.01175**2), dtype=np.longdouble))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    parser.add_argument("--time", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    restart = args.run / "khout.000001.dat"
    start_time = snapshot_time(restart, COUNT)
    if not start_time < args.time:
        raise ValueError("Replay endpoint must follow accepted restart")
    candidates = [p for p in args.run.glob("khout.*.dat")
                  if abs(snapshot_time(p, COUNT)-args.time) < 1e-10]
    if len(candidates) != 1:
        raise ValueError("Missing unique complete replay endpoint")
    endpoint = candidates[0]
    start, _ = ordered(restart, COUNT)
    final, _ = ordered(endpoint, COUNT)
    s0, s1 = state_metrics(start), state_metrics(final)
    r0x = start["x"].astype(float)-3.
    r0y = start["y"].astype(float)-3.
    radius = np.hypot(r0x, r0y)
    disk = (radius >= .5) & (radius < 2.)
    phase = (args.time-start_time)*(radius*radius+.01175**2)**(-.75)
    x_exact = 3. + r0x*np.cos(phase)-r0y*np.sin(phase)
    y_exact = 3. + r0x*np.sin(phase)+r0y*np.cos(phase)
    orbit_max = float(np.max(np.hypot(final["x"][disk].astype(float)-x_exact[disk],
                                      final["y"][disk].astype(float)-y_exact[disk])))
    e0, e1 = energy(start), energy(final)
    energy_relative = abs(e1-e0)/abs(e0)
    log = (args.run/"log").read_text(errors="replace")
    ledgers = [line for line in log.splitlines() if "[RK4E]" in line]
    if not ledgers:
        raise ValueError("No RK4 energy ledger")
    injection = {}
    for key in ("floor_cum", "de_cum", "sfl_cum"):
        match = re.search(rf"\b{key}=([0-9.eE+-]+)", ledgers[-1])
        if match is None:
            raise ValueError(f"Missing {key} in ledger")
        injection[key] = float(match.group(1))
    checks = {
        "positive_state": bool(s0["positive"] and s1["positive"]),
        "fixed_mass_exact": bool(np.array_equal(start["mass"], final["mass"])),
        "W0_exact": bool(np.all(start["w2"] == 0) and np.all(final["w2"] == 0)),
        "bounded_area_tiling": abs(s1["volume_sum"]-DOMAIN_AREA) < 1e-5,
        "density_and_EOS": max(s1["density_relative_max"],
                               s1["eos_pressure_relative_max"]) < 1e-6,
        "energy": energy_relative < 1e-6,
        "disk_incremental_orbit": orbit_max < 1e-4,
        "zero_injection": all(value == 0.0 for value in injection.values()),
        "no_geometry_abort": "[GPU_CLIP_CAPACITY]" not in log and
                             "[GPU_TESS_INVALID]" not in log,
        "no_wall_crossing": "[FIXED_POLYGON_CROSSING]" not in log,
        "fixed_polygon_enabled": "[FIXED_POLYGON] enabled" in log,
    }
    checks = {key: bool(value) for key, value in checks.items()}
    receipt = {
        "run": str(args.run.resolve()), "restart": str(restart.resolve()),
        "restart_sha256": digest(restart), "endpoint": str(endpoint.resolve()),
        "endpoint_sha256": digest(endpoint), "start_time": start_time,
        "endpoint_time": args.time, "initial_state": s0, "final_state": s1,
        "energy_relative_change_from_restart": energy_relative,
        "disk_incremental_orbit_error_max": orbit_max,
        "rk4_injection": injection, "checks": checks,
        "scope": "A short accepted-state replay of the revised tessellation; no t>6 stability claim",
    }
    rendered = json.dumps(receipt, indent=2) + "\n"
    with args.output.open("x") as stream:
        stream.write(rendered)
    print(rendered)
    if not all(checks.values()):
        raise SystemExit("Bounded short replay gate failed")


if __name__ == "__main__":
    main()

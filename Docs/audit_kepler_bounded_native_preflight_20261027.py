"""Gate the first native CUDA W=0 fixed-polygon Kepler step."""

import argparse
import json
from pathlib import Path
import re
import struct

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_state_consistency_20261027 import state_metrics
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N = 55351
GAMMA_MINUS_ONE = float(np.float32(1.666667)) - 1.0
POLYGON_AREA = 0.5 * 1024 * 2.75**2 * np.sin(2 * np.pi / 1024)


def profile(radius):
    rho = np.empty_like(radius)
    inner = radius < .5
    outer = radius > 2.
    rho[inner] = .01 + (radius[inner] / .5)**3
    rho[~inner & ~outer] = 1.01
    rho[outer] = .01 + (1. + (radius[outer] - 2.) / .1)**-3
    return rho


def energy(p):
    x = p["x"].astype(float)-3.
    y = p["y"].astype(float)-3.
    m = p["mass"].astype(float)
    vx, vy = p["vx"].astype(float), p["vy"].astype(float)
    return float(np.sum(.5*m*(vx*vx+vy*vy) + p["ie"].astype(float)
                        - m / np.sqrt(x*x+y*y+.01175**2), dtype=np.longdouble))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    run = args.run
    initial, endpoint = run/"khout.000000.dat", run/"khout.000001.dat"
    t0, t1 = snapshot_time(initial, N), snapshot_time(endpoint, N)
    if abs(t0) > 1e-12 or abs(t1-.001) > 1e-12:
        raise ValueError("Missing exact complete native preflight states")
    p0, _ = ordered(initial, N)
    p1, _ = ordered(endpoint, N)
    raw = (run/"areas.bin").read_bytes()
    if len(raw) != 20+8*N or struct.unpack_from("<idd", raw) != (N,6.,6.):
        raise ValueError("Invalid bounded-area sidecar")
    areas = np.frombuffer(raw, dtype="<f8", count=N, offset=20)
    r0 = np.hypot(p0["x"].astype(float)-3., p0["y"].astype(float)-3.)
    disk = (r0>=.5)&(r0<2.)
    m0 = p0["mass"].astype(float)
    initial_state, final_state = state_metrics(p0), state_metrics(p1)
    expected_mass = profile(r0)*areas
    mass_rel = float(np.max(np.abs(m0/expected_mass-1.)))
    initial_area_rel = float(np.max(np.abs(p0["volume"].astype(float)/areas-1.)))
    x0, y0 = p0["x"].astype(float)-3., p0["y"].astype(float)-3.
    phase = t1*(r0*r0+.01175**2)**(-.75)
    x_exact = 3.+x0*np.cos(phase)-y0*np.sin(phase)
    y_exact = 3.+x0*np.sin(phase)+y0*np.cos(phase)
    disk_orbit_max = float(np.max(np.hypot(
        p1["x"][disk].astype(float)-x_exact[disk],
        p1["y"][disk].astype(float)-y_exact[disk])))
    log = (run/"log").read_text(errors="replace")
    ledger_lines = [line for line in log.splitlines() if "[RK4E]" in line]
    if not ledger_lines:
        raise ValueError("Missing native RK4 energy ledger")
    last = ledger_lines[-1]
    ledger = {}
    for key in ("floor_cum", "de_cum", "sfl_cum"):
        match = re.search(rf"\b{key}=([0-9.eE+-]+)", last)
        if not match:
            raise ValueError(f"Missing {key}")
        ledger[key] = float(match.group(1))
    e0, e1 = energy(p0), energy(p1)
    e_rel = abs(e1-e0)/abs(e0)
    checks = {
        "sidecar_tiles_fixed_polygon": abs(float(np.sum(areas,dtype=np.longdouble))-POLYGON_AREA)<1e-8,
        "initial_native_areas_match_sidecar": initial_area_rel<1e-6,
        "initial_native_area_tiles_polygon": abs(initial_state["volume_sum"]-POLYGON_AREA)<1e-6,
        "evolved_native_area_tiles_polygon": abs(final_state["volume_sum"]-POLYGON_AREA)<1e-5,
        "initial_mass_matches_profile_and_area": mass_rel<1e-6,
        "fixed_mass_exact": bool(np.array_equal(p0["mass"],p1["mass"])),
        "W0_exact": bool(np.all(p0["w2"]==0)&np.all(p1["w2"]==0)),
        "positive_state": bool(initial_state["positive"] and final_state["positive"]),
        "rho_and_EOS": max(initial_state["density_relative_max"],
                           initial_state["eos_pressure_relative_max"],
                           final_state["density_relative_max"],
                           final_state["eos_pressure_relative_max"])<1e-6,
        "disk_orbit": disk_orbit_max<1e-5,
        "energy": e_rel<1e-4,
        "zero_injection": all(value==0.0 for value in ledger.values()),
        "no_wall_crossing": "[FIXED_POLYGON_CROSSING]" not in log,
        "fixed_polygon_enabled": "[FIXED_POLYGON] enabled" in log,
    }
    checks = {name: bool(value) for name, value in checks.items()}
    receipt = {
        "run": str(run.resolve()), "initial_sha256": digest(initial),
        "endpoint_sha256": digest(endpoint), "area_sha256": digest(run/"areas.bin"),
        "n": N, "t0": t0, "t1": t1, "polygon_area": POLYGON_AREA,
        "sidecar_area_sum": float(np.sum(areas,dtype=np.longdouble)),
        "initial_state": initial_state, "final_state": final_state,
        "initial_area_relative_max": initial_area_rel,
        "initial_mass_relative_max": mass_rel,
        "initial_disk_sites": int(np.sum(disk)),
        "disk_orbit_position_error_max": disk_orbit_max,
        "energy_relative_change": e_rel, "ledger": ledger,
        "checks": checks,
        "scope": "First one-rank native CUDA fixed-polygon W0 step only; no multi-rank or long-time inference",
    }
    with args.output.open("x") as stream:
        json.dump(receipt,stream,indent=2)
        stream.write("\n")
    print(json.dumps(receipt,indent=2))
    if not all(checks.values()):
        raise SystemExit("Bounded native preflight gate failed")


if __name__ == "__main__":
    main()

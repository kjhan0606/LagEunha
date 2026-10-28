"""Independent native sparse Kepler state/thermodynamic consistency gate.

This complements the time-specific orbital and conservation auditors.  It
does not infer long-time stability from a positive snapshot.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import IC, digest, snapshot_time


GAMMA_MINUS_ONE = float(np.float32(1.666667)) - 1.


def state_metrics(p):
    mass = p["mass"].astype(float)
    volume = p["volume"].astype(float)
    internal = p["ie"].astype(float)
    density = p["den"].astype(float)
    pressure = p["pressure"].astype(float)
    positive = bool(np.all(mass > 0) and np.all(volume > 0) and
                    np.all(internal > 0) and np.all(density > 0) and
                    np.all(pressure > 0))
    if not positive:
        return {"positive": False, "density_relative_max": None,
                "eos_pressure_relative_max": None,
                "volume_sum": float(np.sum(volume, dtype=np.longdouble))}
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        density_reference = mass / volume
        pressure_reference = GAMMA_MINUS_ONE * internal / volume
        rho_error = np.abs(density / density_reference - 1.)
        pressure_error = np.abs(pressure / pressure_reference - 1.)
    return {
        "positive": True,
        "density_relative_max": float(np.max(rho_error)),
        "eos_pressure_relative_max": float(np.max(pressure_error)),
        "volume_sum": float(np.sum(volume, dtype=np.longdouble)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("initial", type=Path)
    ap.add_argument("endpoint", type=Path)
    ap.add_argument("--expected-count", type=int, choices=tuple(IC), required=True)
    ap.add_argument("--time", type=float, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    n = args.expected_count
    if not np.isfinite(args.time) or args.time <= 0:
        raise ValueError("Expected finite positive endpoint time")
    t0 = snapshot_time(args.initial, n)
    t1 = snapshot_time(args.endpoint, n)
    if abs(t0) > 1e-12 or abs(t1 - args.time) > 1e-10:
        raise ValueError("Native snapshot times do not match requested gate")
    p0, meta0 = ordered(args.initial, n)
    p1, meta1 = ordered(args.endpoint, n)
    if abs(float(meta0["t"]) - t0) > 1e-12 or abs(float(meta1["t"]) - t1) > 1e-12:
        raise ValueError("Snapshot header/footer time disagreement")
    fields = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den",
              "pressure", "w2")
    finite = bool(all(np.isfinite(p[field]).all() for p in (p0, p1)
                      for field in fields))
    a, b = state_metrics(p0), state_metrics(p1)
    checks = {
        "finite": finite,
        "positive_state": a["positive"] and b["positive"],
        "fixed_mass_exact": bool(np.array_equal(p0["mass"], p1["mass"])),
        "W0_exact": bool(np.all(p0["w2"] == 0) and np.all(p1["w2"] == 0)),
        "initial_area_tiles_box": bool(np.isfinite(a["volume_sum"]) and
                                       abs(a["volume_sum"] - 36.) < 1e-6),
        "endpoint_area_tiles_box": bool(np.isfinite(b["volume_sum"]) and
                                        abs(b["volume_sum"] - 36.) < 1e-5),
        "rho_equals_mass_over_volume": bool(a["density_relative_max"] is not None and
                                             b["density_relative_max"] is not None and
                                             max(a["density_relative_max"],
                                                 b["density_relative_max"]) < 1e-6),
        "pressure_matches_ideal_gas_EOS": bool(a["eos_pressure_relative_max"] is not None and
                                                 b["eos_pressure_relative_max"] is not None and
                                                 max(a["eos_pressure_relative_max"],
                                                     b["eos_pressure_relative_max"]) < 1e-6),
    }
    receipt = {
        "initial": str(args.initial.resolve()),
        "initial_sha256": digest(args.initial),
        "endpoint": str(args.endpoint.resolve()),
        "endpoint_sha256": digest(args.endpoint),
        "n": n, "time": t1, "gamma_minus_one": GAMMA_MINUS_ONE,
        "initial_state": a, "endpoint_state": b, "checks": checks,
        "scope": "Native W0 sparse state consistency; requires separate orbital, energy, entropy and boundary audits",
    }
    with args.output.open("x") as stream:
        json.dump(receipt, stream, indent=2)
        stream.write("\n")
    print(json.dumps(receipt, indent=2))
    if not all(checks.values()):
        raise SystemExit("Native sparse state consistency gate failed")


if __name__ == "__main__":
    main()

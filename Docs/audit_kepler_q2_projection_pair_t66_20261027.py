"""Audit same-GPU projected/raw moving weights through native t=6.6."""

import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_bounded_native_evolution_20261027 import energy
from audit_kepler_q2_raw_rate_pair_t6501_20261027 import STAGE, metrics
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
RUNS = {
    "projected": ROOT / "gpu_energy_neutral_projected_v11b_t65_to66_ada_v1",
    "raw": ROOT / "gpu_energy_neutral_raw_rate_v11b_t65_to66_ada_v1",
}
OUTPUT = Path(__file__).resolve().parent / "KEPLER_Q2_PROJECTION_NATIVE_PAIR_T66_2026-10-27.json"
TIME = re.compile(r"Time=\s*([\deE+.-]+)\s+step=\s*(\d+)\s+dt=\s*([\deE+.-]+)")


def work_ledger(run, expected_time=6.6):
    stage = []
    steps = []
    for line in (run / "log").read_text(errors="replace").splitlines():
        match = STAGE.fullmatch(line)
        if match:
            stage.append((int(match.group(1)), float(match.group(4)),
                          float(match.group(5))))
        clock = TIME.fullmatch(line)
        if clock:
            if len(stage) != 4 or [item[0] for item in stage] != list(range(4)):
                raise ValueError(f"Incomplete native RK4 weight stages before {line}")
            dt = float(clock.group(3))
            weights = np.array((1., 2., 2., 1.))/6.
            steps.append({"time": float(clock.group(1)),
                          "step": int(clock.group(2)), "dt": dt,
                          "integrated_weight_work": float(dt*np.dot(
                              weights, [item[2] for item in stage])),
                          "max_abs_alpha": max(abs(item[1]) for item in stage),
                          "max_abs_work_rate": max(abs(item[2]) for item in stage)})
            stage.clear()
    if stage or not steps or abs(steps[-1]["time"]-expected_time) > 1e-10:
        raise ValueError(f"Incomplete native t={expected_time} work ledger")
    return {"steps": len(steps),
            "integrated_weight_work": float(sum(row["integrated_weight_work"]
                                                for row in steps)),
            "max_abs_alpha": max(row["max_abs_alpha"] for row in steps),
            "max_abs_work_rate": max(row["max_abs_work_rate"] for row in steps),
            "min_dt": min(row["dt"] for row in steps),
            "max_dt": max(row["dt"] for row in steps),
            "first_step_work": steps[0]["integrated_weight_work"],
            "last_step_work": steps[-1]["integrated_weight_work"]}


def main():
    receipts = {name: json.loads((run / "audit_t66.json").read_text())
                for name, run in RUNS.items()}
    a, b = receipts["projected"], receipts["raw"]
    if (a["time"] != 6.6 or b["time"] != 6.6
            or a["parent_time"] != 6.5 or b["parent_time"] != 6.5
            or a["initial_sha256"] != b["initial_sha256"]
            or a["parent_sha256"] != b["parent_sha256"]
            or a["binary_sha256"] != b["binary_sha256"]
            or not all(a["checks"].values()) or not all(b["checks"].values())):
        raise ValueError("Native t=6.6 pair has invalid lineage or physics gate")
    for receipt in receipts.values():
        for key, sha_key in (("initial", "initial_sha256"),
                            ("parent", "parent_sha256"),
                            ("endpoint", "endpoint_sha256")):
            if digest(Path(receipt[key])) != receipt[sha_key]:
                raise ValueError(f"Changed native {key} snapshot")
    n = a["n"]
    initial, _ = ordered(Path(a["initial"]), n)
    start, _ = ordered(Path(a["parent"]), n)
    endpoints = {name: ordered(Path(receipt["endpoint"]), n)[0]
                 for name, receipt in receipts.items()}
    r0 = np.hypot(initial["x"].astype(float)-3.,
                  initial["y"].astype(float)-3.)
    ledgers = {name: work_ledger(run) for name, run in RUNS.items()}
    if ledgers["raw"]["max_abs_alpha"] != 0.:
        raise ValueError("Raw rate was unexpectedly projected")
    e0 = energy(start)
    result = {"time_start": 6.5, "time_end": 6.6,
              "binary_sha256": a["binary_sha256"],
              "shared_parent_sha256": a["parent_sha256"],
              "receipts": {name: str(run / "audit_t66.json")
                           for name, run in RUNS.items()},
              "endpoint_sha256": {name: receipt["endpoint_sha256"]
                                  for name, receipt in receipts.items()},
              "ledgers": ledgers,
              "energy_parent": e0,
              "energy_change": {name: energy(endpoint)-e0
                                for name, endpoint in endpoints.items()},
              "energy_raw_minus_projected": (energy(endpoints["raw"])
                                             - energy(endpoints["projected"])),
              "bins": {name: metrics(initial, start, endpoint, r0)
                       for name, endpoint in endpoints.items()},
              "scope": "Same Ada GPU, v11b binary and accepted t=6.5 state. "
                       "First-ring pressure and torque "
                       "may be inherited from the parent; this is not a "
                       "t=0 causal prevention test."}
    with OUTPUT.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
    for key in ("0.5:0.52", "0.52:0.55", "0.5:2"):
        p, raw = (result["bins"][mode][key] for mode in ("projected", "raw"))
        print(key, "volume_rms projected/raw", p["volume_step_rms"],
              raw["volume_step_rms"], "pressure_error projected/raw",
              p["pressure_error_mean_abs"], raw["pressure_error_mean_abs"],
              "ell_increment projected/raw", p["ell_step_mean_abs_relative_initial"],
              raw["ell_step_mean_abs_relative_initial"])
    print("weight_work", ledgers["projected"]["integrated_weight_work"],
          ledgers["raw"]["integrated_weight_work"],
          "energy_raw_minus_projected", result["energy_raw_minus_projected"])


if __name__ == "__main__":
    main()

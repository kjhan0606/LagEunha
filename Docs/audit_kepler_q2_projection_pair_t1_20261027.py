"""Audit native raw/projected early pressure growth from t=.825 to t=1."""

import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_native_evolution_20261027 import energy
from audit_kepler_q2_projection_pair_t66_20261027 import work_ledger
from audit_kepler_q2_raw_rate_pair_t6501_20261027 import metrics
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
RUNS = {
    "projected": ROOT / "gpu_energy_neutral_projected_v11b_t0825_to1_ada_v1",
    "raw": ROOT / "gpu_energy_neutral_raw_rate_v11b_t0825_to1_ada_v1",
}
OUTPUT = Path(__file__).resolve().parent / "KEPLER_Q2_PROJECTION_NATIVE_EARLY_PAIR_T1_2026-10-27.json"


def main():
    audits = {mode: json.loads((run / "audit_t1.json").read_text())
              for mode, run in RUNS.items()}
    a, b = audits["projected"], audits["raw"]
    if (a["time"] != 1. or b["time"] != 1.
            or a["parent_time"] != .825 or b["parent_time"] != .825
            or a["initial_sha256"] != b["initial_sha256"]
            or a["parent_sha256"] != b["parent_sha256"]
            or a["binary_sha256"] != b["binary_sha256"]
            or not all(a["checks"].values()) or not all(b["checks"].values())):
        raise ValueError("Early native pair lineage or physics gate failed")
    for receipt in audits.values():
        for key, sha_key in (("initial", "initial_sha256"),
                            ("parent", "parent_sha256"),
                            ("endpoint", "endpoint_sha256")):
            if digest(Path(receipt[key])) != receipt[sha_key]:
                raise ValueError(f"Changed native {key} snapshot")
    initial, _ = ordered(Path(a["initial"]), a["n"])
    parent, _ = ordered(Path(a["parent"]), a["n"])
    endpoints = {mode: ordered(Path(receipt["endpoint"]), a["n"])[0]
                 for mode, receipt in audits.items()}
    r0 = np.hypot(initial["x"].astype(float)-3.,
                  initial["y"].astype(float)-3.)
    ledgers = {mode: work_ledger(run, 1.) for mode, run in RUNS.items()}
    if ledgers["raw"]["max_abs_alpha"] != 0.:
        raise ValueError("Raw rate was unexpectedly projected")
    result = {"time_start": .825, "time_end": 1.,
              "initial_sha256": a["initial_sha256"],
              "shared_parent_sha256": a["parent_sha256"],
              "binary_sha256": a["binary_sha256"],
              "receipts": {mode: str(run / "audit_t1.json")
                           for mode, run in RUNS.items()},
              "endpoint_sha256": {mode: receipt["endpoint_sha256"]
                                  for mode, receipt in audits.items()},
              "ledgers": ledgers,
              "energy_parent": energy(parent),
              "energy_change": {mode: energy(endpoint)-energy(parent)
                                for mode, endpoint in endpoints.items()},
              "energy_raw_minus_projected": (energy(endpoints["raw"])
                                             -energy(endpoints["projected"])),
              "bins": {mode: metrics(initial, parent, endpoint, r0)
                       for mode, endpoint in endpoints.items()},
              "scope": "Same accepted t=.825 checkpoint, Ada GPU and v11b "
                       "binary. Earlier t<.825 pressure noise is inherited; "
                       "this is a prevention test only for the following .175 time."}
    with OUTPUT.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
    for key in ("0.5:0.52", "0.5:2"):
        p, raw = (result["bins"][mode][key] for mode in ("projected", "raw"))
        print(key, "volume_rms projected/raw", p["volume_step_rms"],
              raw["volume_step_rms"], "pressure projected/raw",
              p["pressure_error_mean_abs"], raw["pressure_error_mean_abs"],
              "ell_increment projected/raw", p["ell_step_mean_abs_relative_initial"],
              raw["ell_step_mean_abs_relative_initial"])
    print("raw_weight_work", ledgers["raw"]["integrated_weight_work"],
          "energy_raw_minus_projected", result["energy_raw_minus_projected"])


if __name__ == "__main__":
    main()

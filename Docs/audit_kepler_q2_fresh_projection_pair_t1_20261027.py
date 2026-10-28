"""Audit native same-t0 raw/projected q2 trajectories at exact t=1."""

import json
from pathlib import Path

import numpy as np

from audit_kepler_q2_fresh_projection_pair_t001_20261027 import ring_metrics
from audit_kepler_q2_projection_pair_t66_20261027 import work_ledger
from audit_kepler_q2_raw_rate_pair_t6501_20261027 import metrics
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PILOT_PAIR = DOCS / "KEPLER_Q2_FRESH_PROJECTION_PAIR_T001_2026-10-27.json"
OUTPUT = DOCS / "KEPLER_Q2_FRESH_PROJECTION_PAIR_T1_2026-10-27.json"
MODES = ("projected", "raw")


def main():
    pilot = json.loads(PILOT_PAIR.read_text())
    result = {"time_start": 0., "time_end": 1., "binary_sha256": None,
              "initial_sha256": None, "modes": {},
              "scope": "Independent native raw/projected evolution from the same fresh t=0 "
                       "jitter IC. This tests initial noise generation, not t>6 stability."}
    for mode in MODES:
        run = ROOT / f"gpu_energy_neutral_v11b_fresh_{mode}_t001_to1_a10_v1"
        manifest = json.loads((run/"stage_manifest.json").read_text())
        audit = json.loads((run/"audit_t1.json").read_text())
        if (manifest["mode"] != mode or audit["time"] != 1.
                or audit["parent_time"] != .001
                or audit["binary_sha256"] != manifest["binary_sha256"]
                or audit["initial_sha256"] != manifest["initial_sha256"]
                or audit["parent_sha256"] != manifest["accepted_parent_endpoint_sha256"]
                or audit["parent_sha256"] != pilot["modes"][mode]["endpoint_sha256"]
                or len(audit["checks"]) != 21 or not all(audit["checks"].values())):
            raise ValueError(f"Unaccepted fresh {mode} t=1 native trajectory")
        for key, sha_key in (("initial","initial_sha256"),
                            ("parent","parent_sha256"),
                            ("endpoint","endpoint_sha256")):
            if digest(Path(audit[key])) != audit[sha_key]:
                raise ValueError(f"Changed fresh {mode} {key} snapshot")
        if result["initial_sha256"] is None:
            result["initial_sha256"] = audit["initial_sha256"]
            result["binary_sha256"] = audit["binary_sha256"]
        elif (audit["initial_sha256"] != result["initial_sha256"] or
              audit["binary_sha256"] != result["binary_sha256"]):
            raise ValueError("Fresh t=1 pair differs in binary or accepted t=0 reference")
        initial,_ = ordered(Path(audit["initial"]), audit["n"])
        parent,_ = ordered(Path(audit["parent"]), audit["n"])
        endpoint,_ = ordered(Path(audit["endpoint"]), audit["n"])
        x0,y0=initial["x"].astype(float)-3,initial["y"].astype(float)-3
        r0=np.hypot(x0,y0)
        disk=(r0>=.5)&(r0<2)
        ring=(r0>=.5)&(r0<.52)
        ring_ids=np.flatnonzero(ring)
        order=ring_ids[np.argsort(np.arctan2(y0[ring_ids],x0[ring_ids]))]
        ledger = work_ledger(run,1.)
        if mode == "raw" and ledger["max_abs_alpha"] != 0.:
            raise ValueError("Fresh raw branch unexpectedly applied projection")
        pilot_work = pilot["modes"][mode]["gpu_stage_pressure_weight_works"]
        pilot_integral=.001*float(np.dot([1.,2.,2.,1.],pilot_work)/6.)
        result["modes"][mode]={
            "run":str(run),"receipt":str(run/"audit_t1.json"),
            "endpoint_sha256":audit["endpoint_sha256"],
            "parent_sha256":audit["parent_sha256"],
            "checks_passed":len(audit["checks"]),
            "pilot_integrated_weight_work":pilot_integral,
            "continuation_ledger":ledger,
            "total_integrated_weight_work":pilot_integral+ledger["integrated_weight_work"],
            "radial_bins":metrics(initial,parent,endpoint,r0),
            **ring_metrics(initial,endpoint,ring,order,disk),
        }
    p,r=(result["modes"][m] for m in MODES)
    keys=("first_ring_volume_change_rms","first_ring_density_change_rms",
          "first_ring_pressure_change_mean_abs","first_ring_angular_pressure_jump_rms",
          "first_ring_orbital_ell_change_mean_abs")
    result["ratios_projected_over_raw"]={key:p[key]/r[key] for key in keys if r[key]!=0}
    with OUTPUT.open("x") as output:
        json.dump(result,output,indent=2,allow_nan=False)
        output.write("\n")
    print(json.dumps({"output":str(OUTPUT),"ratios":result["ratios_projected_over_raw"],
                      "raw_total_weight_work":r["total_integrated_weight_work"]},indent=2))


if __name__ == "__main__":
    main()

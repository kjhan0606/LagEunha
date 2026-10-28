"""Compare accepted q2 raw-from-.825 trajectory with projected t=6.5 baseline."""

import json
from pathlib import Path

import numpy as np

from audit_kepler_q2_fresh_projection_pair_t001_20261027 import ring_metrics
from audit_kepler_q2_projection_pair_t66_20261027 import work_ledger
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DOCS = Path(__file__).resolve().parent
BASE = ROOT / "gpu_energy_neutral_capacity64_v5_recover_to65_a10_v1"
RAW_EARLY = ROOT / "gpu_energy_neutral_raw_rate_v11b_t0825_to1_ada_v1"
RAW_LONG = ROOT / "gpu_energy_neutral_raw_rate_v11b_t1_to65_ada_v1"
OUTPUT = DOCS / "KEPLER_Q2_RAW_FROM0825_VS_PROJECTED_T65_2026-10-27.json"


def load_receipt(run, name, target):
    d=json.loads((run/name).read_text())
    if (d["time"]!=target or len(d["checks"])!=21
            or not all(d["checks"].values())):
        raise ValueError(f"Unaccepted {run/name}")
    for key,sha in (("initial","initial_sha256"),
                    ("parent","parent_sha256"),
                    ("endpoint","endpoint_sha256")):
        if digest(Path(d[key]))!=d[sha]:
            raise ValueError(f"Changed {key}: {run/name}")
    return d


def main():
    base=load_receipt(BASE,"audit_t65.json",6.5)
    early=load_receipt(RAW_EARLY,"audit_t1.json",1.)
    raw=load_receipt(RAW_LONG,"audit_t65.json",6.5)
    if (base["parent_time"]!=.825 or early["parent_time"]!=.825
            or raw["parent_time"]!=1.
            or base["initial_sha256"]!=early["initial_sha256"]!=raw["initial_sha256"]
            or base["parent_sha256"]!=early["parent_sha256"]
            or raw["parent_sha256"]!=early["endpoint_sha256"]
            or raw["binary_sha256"]!=early["binary_sha256"]):
        raise ValueError("Raw/projected long branches do not share accepted lineage")
    initial,_=ordered(Path(base["initial"]),base["n"])
    x0,y0=initial["x"].astype(float)-3.,initial["y"].astype(float)-3.
    r0=np.hypot(x0,y0)
    ring=(r0>=.5)&(r0<.52)
    disk=(r0>=.5)&(r0<2.)
    ring_ids=np.flatnonzero(ring)
    order=ring_ids[np.argsort(np.arctan2(y0[ring_ids],x0[ring_ids]))]
    metrics={}
    for mode,receipt in (("projected",base),("raw",raw)):
        end,_=ordered(Path(receipt["endpoint"]),receipt["n"])
        volume=end["volume"][ring].astype(float)
        ie=end["ie"][ring].astype(float)
        entropy=np.log(ie/initial["ie"][ring].astype(float))
        entropy+=(float(np.float32(1.666667))-1.)*np.log(
            volume/initial["volume"][ring].astype(float))
        metrics[mode]={"receipt":str((BASE if mode=="projected" else RAW_LONG)
                                      /"audit_t65.json"),
                       "endpoint_sha256":receipt["endpoint_sha256"],
                       "binary_sha256":receipt["binary_sha256"],
                       "first_ring_entropy_log_mean_abs":float(np.mean(np.abs(entropy))),
                       **ring_metrics(initial,end,ring,order,disk)}
    raw_early_work=work_ledger(RAW_EARLY,1.)
    raw_long_work=work_ledger(RAW_LONG,6.5)
    if raw_early_work["max_abs_alpha"]!=0 or raw_long_work["max_abs_alpha"]!=0:
        raise ValueError("Raw long branch unexpectedly applied projection")
    p,r=metrics["projected"],metrics["raw"]
    keys=("first_ring_volume_change_rms","first_ring_density_change_rms",
          "first_ring_pressure_change_mean_abs","first_ring_angular_pressure_jump_rms",
          "first_ring_orbital_ell_change_mean_abs","first_ring_entropy_log_mean_abs")
    result={"shared_initial_sha256":base["initial_sha256"],
            "shared_t0825_sha256":base["parent_sha256"],
            "raw_t1_sha256":early["endpoint_sha256"],
            "time":6.5,"modes":metrics,
            "ratios_projected_over_raw":{k:p[k]/r[k] for k in keys if r[k]!=0},
            "raw_pressure_weight_work_t0825_to65":
                raw_early_work["integrated_weight_work"]+raw_long_work["integrated_weight_work"],
            "raw_work_ledgers":{"early":raw_early_work,"long":raw_long_work},
            "scope":"Both branches share exact t=0 and t=.825 states. Raw uses diagnostic "
                    "v11b from .825; projected uses accepted v5 from .825. Default v11b "
                    "was independently verified bitwise against v5 for one step, but the "
                    "full t=6.5 pair is not the same executable."}
    with OUTPUT.open("x") as output:
        json.dump(result,output,indent=2,allow_nan=False)
        output.write("\n")
    print(json.dumps({"output":str(OUTPUT),"ratios":result["ratios_projected_over_raw"],
                      "raw_pressure_weight_work":result["raw_pressure_weight_work_t0825_to65"]},indent=2))


if __name__=="__main__":
    main()

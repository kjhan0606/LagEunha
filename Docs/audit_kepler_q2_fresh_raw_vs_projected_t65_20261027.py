"""Compare accepted fresh-from-t0 raw q2 trajectory with projected t=6.5."""

import json
from pathlib import Path

import numpy as np

from audit_kepler_q2_fresh_projection_pair_t001_20261027 import ring_metrics
from audit_kepler_q2_projection_pair_t66_20261027 import work_ledger
from audit_kepler_q2_raw_from0825_vs_projected_t65_20261027 import load_receipt
from audit_kepler_sparse_area_t0825_20261027 import ordered


ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DOCS=Path(__file__).resolve().parent
BASE=ROOT/"gpu_energy_neutral_capacity64_v5_recover_to65_a10_v1"
RAW_PILOT=ROOT/"gpu_energy_neutral_v11b_fresh_raw_t001_a10_v1"
RAW_T1=ROOT/"gpu_energy_neutral_v11b_fresh_raw_t001_to1_a10_v1"
RAW_T65=ROOT/"gpu_energy_neutral_raw_rate_v11b_fresh_t1_to65_ada_v1"
PILOT_PAIR=DOCS/"KEPLER_Q2_FRESH_PROJECTION_PAIR_T001_2026-10-27.json"
OUTPUT=DOCS/"KEPLER_Q2_FRESH_RAW_VS_PROJECTED_T65_2026-10-27.json"


def main():
    base=load_receipt(BASE,"audit_t65.json",6.5)
    raw_t1=load_receipt(RAW_T1,"audit_t1.json",1.)
    raw_t65=load_receipt(RAW_T65,"audit_t65.json",6.5)
    pilot=json.loads(PILOT_PAIR.read_text())
    pilot_raw=json.loads((RAW_PILOT/"audit_t001.json").read_text())
    if (base["initial_sha256"]!=raw_t1["initial_sha256"]!=raw_t65["initial_sha256"]
            or raw_t1["parent_time"]!=.001 or raw_t65["parent_time"]!=1.
            or raw_t1["parent_sha256"]!=pilot_raw["endpoint_sha256"]
            or raw_t65["parent_sha256"]!=raw_t1["endpoint_sha256"]
            or raw_t65["binary_sha256"]!=raw_t1["binary_sha256"]
            or pilot["modes"]["raw"]["endpoint_sha256"]!=pilot_raw["endpoint_sha256"]
            or len(pilot_raw["checks"])!=18 or not all(pilot_raw["checks"].values())):
        raise ValueError("Fresh t=0 raw lineage or accepted pilot differs")
    initial,_=ordered(Path(base["initial"]),base["n"])
    x0,y0=initial["x"].astype(float)-3.,initial["y"].astype(float)-3.
    r0=np.hypot(x0,y0)
    ring=(r0>=.5)&(r0<.52)
    disk=(r0>=.5)&(r0<2.)
    ids=np.flatnonzero(ring)
    order=ids[np.argsort(np.arctan2(y0[ids],x0[ids]))]
    rows={}
    for mode,receipt,run in (("projected",base,BASE),("raw",raw_t65,RAW_T65)):
        end,_=ordered(Path(receipt["endpoint"]),receipt["n"])
        entropy=np.log(end["ie"][ring].astype(float)/initial["ie"][ring].astype(float))
        entropy+=(float(np.float32(1.666667))-1.)*np.log(
            end["volume"][ring].astype(float)/initial["volume"][ring].astype(float))
        rows[mode]={"receipt":str(run/"audit_t65.json"),
                    "binary_sha256":receipt["binary_sha256"],
                    "endpoint_sha256":receipt["endpoint_sha256"],
                    "first_ring_entropy_log_mean_abs":float(np.mean(np.abs(entropy))),
                    **ring_metrics(initial,end,ring,order,disk)}
    pilot_works=pilot["modes"]["raw"]["gpu_stage_pressure_weight_works"]
    pilot_work=.001*float(np.dot([1.,2.,2.,1.],pilot_works)/6.)
    work_t1=work_ledger(RAW_T1,1.)
    work_t65=work_ledger(RAW_T65,6.5)
    if work_t1["max_abs_alpha"]!=0. or work_t65["max_abs_alpha"]!=0.:
        raise ValueError("Fresh raw trajectory unexpectedly projected weights")
    p,r=rows["projected"],rows["raw"]
    keys=("first_ring_volume_change_rms","first_ring_density_change_rms",
          "first_ring_pressure_change_mean_abs","first_ring_angular_pressure_jump_rms",
          "first_ring_orbital_ell_change_mean_abs","first_ring_entropy_log_mean_abs")
    result={"time":6.5,"shared_t0_sha256":base["initial_sha256"],
            "raw_t001_sha256":pilot_raw["endpoint_sha256"],
            "raw_t1_sha256":raw_t1["endpoint_sha256"],"modes":rows,
            "ratios_projected_over_raw":{k:p[k]/r[k] for k in keys if r[k]!=0},
            "raw_integrated_weight_work_t0_to65":
                pilot_work+work_t1["integrated_weight_work"]+work_t65["integrated_weight_work"],
            "raw_work_ledgers":{"pilot":pilot_work,"to1":work_t1,"to65":work_t65},
            "scope":"Both branches share exact accepted t=0 physical IC, but distinct "
                    "v5/v11b binaries. The raw branch has no projected prehistory. "
                    "Native t=6.5 gates and checkpoint hashes are required before comparison."}
    with OUTPUT.open("x") as output:
        json.dump(result,output,indent=2,allow_nan=False)
        output.write("\n")
    print(json.dumps({"output":str(OUTPUT),"ratios":result["ratios_projected_over_raw"],
                      "raw_weight_work":result["raw_integrated_weight_work_t0_to65"]},indent=2))


if __name__=="__main__":
    main()

"""Measure first-stage projection conditioning for the fresh q2 Kepler IC."""

import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DOCS = Path(__file__).resolve().parent
PAIR = DOCS / "KEPLER_Q2_FRESH_PROJECTION_PAIR_T001_2026-10-27.json"
OUTPUT = DOCS / "KEPLER_Q2_FRESH_PROJECTION_NEARNULL_T001_2026-10-27.json"
STAGE = re.compile(r"^\[GPU_VOLW2\] stage=(\d+) active=(\d+) faces=(\d+) "
                   r"alpha=([\deE+.-]+) work=([\deE+.-]+) elapsed=([\deE+.-]+) s$",re.M)


def main():
    pair=json.loads(PAIR.read_text())
    if pair["time"]!=.001 or len(pair["initial_state_fields_exact"])!=16:
        raise ValueError("Fresh causal pair is not accepted")
    rows={}
    for mode in ("projected","raw"):
        run=ROOT/f"gpu_energy_neutral_v11b_fresh_{mode}_t001_a10_v1"
        audit=json.loads((run/"audit_t001.json").read_text())
        if (len(audit["checks"])!=18 or not all(audit["checks"].values())
                or audit["binary_sha256"]!=pair["modes"][mode]["binary_sha256"]
                or audit["endpoint_sha256"]!=pair["modes"][mode]["endpoint_sha256"]
                or digest(run/"khout.000000.dat")!=audit["initial_sha256"]):
            raise ValueError(f"Unaccepted {mode} native stage")
        log=(run/"log").read_text(errors="replace")
        stages=[{"stage":int(x[0]),"active":int(x[1]),"faces":int(x[2]),
                 "alpha":float(x[3]),"work_rate":float(x[4])}
                for x in STAGE.findall(log)]
        if len(stages)!=4 or [x["stage"] for x in stages]!=list(range(4)):
            raise ValueError(f"Incomplete {mode} weight stages")
        rows[mode]={"run":str(run),"initial_sha256":audit["initial_sha256"],
                    "log_sha256":digest(run/"log"),"stages":stages}
    if rows["projected"]["stages"][0]["faces"]!=rows["raw"]["stages"][0]["faces"]:
        raise ValueError("Different first-stage face counts")
    p0,_=ordered(ROOT/"jitter005_t001_a10_v1/khout.000000.dat",55351)
    pressure=p0["pressure"].astype(float)
    if len(np.unique(pressure))!=1:
        raise ValueError("Accepted t=0 pressure is not uniform")
    gamma_minus_one=float(np.float32(1.666667))-1.
    reconstructed=(p0["ie"].astype(float)/p0["volume"].astype(float)
                   *gamma_minus_one).astype(np.float32)
    stage_p=rows["projected"]["stages"][0]
    stage_r=rows["raw"]["stages"][0]
    if stage_r["alpha"]!=0 or stage_p["alpha"]==0:
        raise ValueError("No first-stage projected/raw distinction")
    inferred_norm2=stage_r["work_rate"]/stage_p["alpha"]
    if not inferred_norm2>0:
        raise ValueError("Inconsistent first-stage h norm estimate")
    result={"initial_reference_sha256":digest(ROOT/"jitter005_t001_a10_v1/khout.000000.dat"),
            "initial_pressure_unique_count":1,"initial_pressure_value":float(pressure[0]),
            "ie_over_volume_pressure_float32_unique_count":int(len(np.unique(reconstructed))),
            "ie_over_volume_pressure_float32_span":float(reconstructed.max()-reconstructed.min()),
            "gamma_minus_one_used":gamma_minus_one,"modes":rows,
            "stage0_inferred_h_dot_h":inferred_norm2,
            "stage0_inferred_h_norm":float(np.sqrt(inferred_norm2)),
            "stage0_inferred_rate_correction_norm":float(abs(stage_p["alpha"])*np.sqrt(inferred_norm2)),
            "inference_scope":"The h norms assume both modes see identical stage-0 GPU geometry, "
                              "pressure gradient and raw rate. Equal t=0 parsed fields and face "
                              "counts support but do not prove that assumption. Reconstructed "
                              "ie/volume pressure levels do not capture actual device-stage pressure."}
    with OUTPUT.open("x") as output:
        json.dump(result,output,indent=2,allow_nan=False)
        output.write("\n")
    print(json.dumps({"output":str(OUTPUT),"h2_inferred":inferred_norm2,
                      "reconstructed_pressure_levels":result["ie_over_volume_pressure_float32_unique_count"]},indent=2))


if __name__=="__main__":
    main()

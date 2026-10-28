"""Check experimental weighted CUDA binary reduces bitwise to accepted W0 gas."""

import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
W0=ROOT/"faceforce_pair_t6501_a10_v1"
NEW=ROOT/"gpu_energy_neutral_v2_w0_reduction_t6501_a10_v1"
OUTPUT=Path(__file__).with_name("KEPLER_Q2_ENERGY_NEUTRAL_V2_W0_REDUCTION_T6501_2026-10-27.json")
FIELDS=("x","y","vx","vy","mass","volume","ie","den","pressure","w2")


def main():
    result={"reference_run":str(W0),"candidate_run":str(NEW),"time":6.501,"cases":{},
            "scope":"Mode-off opt-in weighted CUDA binary versus original accepted W0 binary, same t=6.5 parents. Bitwise physical state reduction only."}
    for label in ("aligned","jitter"):
        old=json.loads((W0/label/"audit_t6501.json").read_text())
        new=json.loads((NEW/label/"audit_t6501.json").read_text())
        if not all(old["checks"].values()) or not all(new["checks"].values()):
            raise ValueError(f"{label}: native one-step physics gate failed")
        if old["restart_sha256"]!=new["restart_sha256"]:
            raise ValueError(f"{label}: restarts differ")
        for audit in (old,new):
            if (digest(Path(audit["endpoint"]))!=audit["endpoint_sha256"] or
                    abs(snapshot_time(Path(audit["endpoint"]),55351)-6.501)>1e-10):
                raise ValueError(f"{label}: endpoint changed")
        p0,_=ordered(Path(old["endpoint"]),55351)
        p1,_=ordered(Path(new["endpoint"]),55351)
        fields={key:{"bitwise_equal":bool(np.array_equal(p0[key],p1[key])),
                     "max_absolute_difference":float(np.max(np.abs(
                         p0[key].astype(float)-p1[key].astype(float))))}
                for key in FIELDS}
        result["cases"][label]={"restart_sha256":old["restart_sha256"],
                                "original_endpoint_sha256":old["endpoint_sha256"],
                                "candidate_endpoint_sha256":new["endpoint_sha256"],
                                "physical_fields":fields,
                                "all_fields_bitwise_equal":bool(all(v["bitwise_equal"]
                                                                  for v in fields.values()))}
    result["checks"]={label:row["all_fields_bitwise_equal"]
                      for label,row in result["cases"].items()}
    with OUTPUT.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(result,indent=2,allow_nan=False))
    if not all(result["checks"].values()):
        raise SystemExit("Weighted binary failed W0 reduction")


if __name__=="__main__": main()

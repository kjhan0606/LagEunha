"""Compare accepted weighted and W0 native density drift over the same RK4 step."""

import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
WEIGHTED=ROOT/"gpu_energy_neutral_weighted_t6501_a10_v2"
W0=ROOT/"faceforce_pair_t6501_a10_v1"
OUTPUT=Path(__file__).with_name("KEPLER_Q2_WEIGHTED_NATIVE_V2_VS_W0_T6501_2026-10-27.json")


def metric(p0,p1,selection):
    mass=p0["mass"][selection].astype(float)
    change=p1["den"][selection].astype(float)/p0["den"][selection].astype(float)-1.
    return {"cells":int(np.sum(selection)),
            "mass_weighted_abs_density_change":float(
                np.sum(np.abs(change)*mass,dtype=np.longdouble)
                /np.sum(mass,dtype=np.longdouble)),
            "abs_density_change_p99":float(np.quantile(np.abs(change),.99))}


def main():
    result={"time_start":6.5,"time_end":6.501,"cases":{},
            "scope":"One native accepted RK4 step from identical t=6.5 states; density drift only, no long-time stability claim."}
    for label in ("aligned","jitter"):
        a=WEIGHTED/label
        b=W0/label
        wa=json.loads((a/"audit_t6501.json").read_text())
        wb=json.loads((b/"audit_t6501.json").read_text())
        if not all(wa["checks"].values()) or not all(wb["checks"].values()):
            raise ValueError(f"{label}: one-step gate did not pass")
        if digest(Path(wa["endpoint"]))!=wa["endpoint_sha256"] or digest(Path(wb["endpoint"]))!=wb["endpoint_sha256"]:
            raise ValueError(f"{label}: endpoint changed")
        if digest(Path(wa["restart"]))!=wa["restart_sha256"] or digest(Path(wb["restart"]))!=wb["restart_sha256"]:
            raise ValueError(f"{label}: restart changed")
        if wa["restart_sha256"]!=wb["restart_sha256"]:
            raise ValueError(f"{label}: compared runs have different physical parents")
        p0,_=ordered(Path(wa["restart"]),55351)
        pw,_=ordered(Path(wa["endpoint"]),55351)
        pv,_=ordered(Path(wb["endpoint"]),55351)
        if (abs(snapshot_time(Path(wa["endpoint"]),55351)-6.501)>1e-10
                or abs(snapshot_time(Path(wb["endpoint"]),55351)-6.501)>1e-10):
            raise ValueError(f"{label}: equal-time comparison failed")
        radius=np.hypot(p0["x"].astype(float)-3.,p0["y"].astype(float)-3.)
        bands={}
        for title,lo,hi in (("inner",.5,.75),("disk",.5,2.)):
            select=(radius>=lo)&(radius<hi)
            w0=metric(p0,pv,select)
            weighted=metric(p0,pw,select)
            bands[title]={"radius_range":[lo,hi],"W0":w0,"weighted":weighted,
                          "W0_over_weighted_mean_abs":(
                              w0["mass_weighted_abs_density_change"]
                              /weighted["mass_weighted_abs_density_change"])}
        result["cases"][label]={
            "parent_sha256":wa["restart_sha256"],
            "W0_endpoint_sha256":wb["endpoint_sha256"],
            "weighted_endpoint_sha256":wa["endpoint_sha256"],
            "weight_range":float(np.ptp(pw["w2"].astype(float))),
            "bands":bands,
        }
    with OUTPUT.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":
    main()

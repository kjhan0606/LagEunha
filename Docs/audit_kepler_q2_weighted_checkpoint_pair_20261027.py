"""Read a complete matching W0/weighted Kepler checkpoint during long runs."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
NEW=ROOT/"gpu_energy_neutral_weighted_t10_a10_v1"
W0={"aligned":ROOT/"to10_near32_a10_v1",
    "jitter":ROOT/"jitter005_to10_a10_v1"}
GAMMA_MINUS_ONE=float(np.float32(1.666667))-1.


def weighted_mean(value,mass):
    return float(np.sum(value*mass,dtype=np.longdouble)/np.sum(mass,dtype=np.longdouble))


def measure(p0,p,t,select):
    mass=p0["mass"][select].astype(float)
    volume=p["volume"][select].astype(float)
    den0=p0["den"][select].astype(float)
    den=p["den"][select].astype(float)
    change=den/den0-1.
    entropy=(np.log(p["ie"][select].astype(float)/p0["ie"][select].astype(float))
             +GAMMA_MINUS_ONE*np.log(volume/p0["volume"][select].astype(float)))
    x0=p0["x"][select].astype(float)-3.
    y0=p0["y"][select].astype(float)-3.
    r0=np.hypot(x0,y0)
    phase=t*(r0*r0+.01175**2)**(-.75)
    xe=x0*np.cos(phase)-y0*np.sin(phase)
    ye=x0*np.sin(phase)+y0*np.cos(phase)
    x=p["x"][select].astype(float)-3.
    y=p["y"][select].astype(float)-3.
    orbit_error=np.hypot(x-xe,y-ye)
    return {"cells":int(select.sum()),
            "density_abs_mass_mean":weighted_mean(np.abs(change),mass),
            "density_abs_p99":float(np.quantile(np.abs(change),.99)),
            "entropy_log_abs_mass_mean":weighted_mean(np.abs(entropy),mass),
            "entropy_log_signed_mass_mean":weighted_mean(entropy,mass),
            "orbit_position_error_abs_mass_mean":weighted_mean(orbit_error,mass)}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--step",type=int,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--weighted-only",action="store_true")
    args=ap.parse_args()
    if args.step<=1 or args.output.exists():
        raise ValueError("Expected a new complete intermediate checkpoint/output")
    filename=f"khout.{args.step:06d}.dat"
    result={"step":args.step,"cases":{},
            "scope":("Read-only weighted bounded q2 intermediate checkpoint; not a terminal physics gate or GIZMO-equivalent IC."
                     if args.weighted_only else
                     "Read-only same-time bounded q2 checkpoint comparison; neither a terminal physics gate nor a GIZMO-equivalent IC.")}
    for label in ("aligned","jitter"):
        new=NEW/label/filename
        tnew=snapshot_time(new,55351)
        p0,_=ordered(NEW/label/"initial.dat",55351)
        pw,_=ordered(new,55351)
        if not np.array_equal(pw["mass"],p0["mass"]):
            raise ValueError(f"{label}: mass changed")
        if not args.weighted_only:
            old=W0[label]/filename
            told=snapshot_time(old,55351)
            if abs(tnew-told)>1e-8:
                raise ValueError(f"{label}: unequal physical times")
            pv,_=ordered(old,55351)
            if not np.array_equal(pv["mass"],p0["mass"]):
                raise ValueError(f"{label}: W0 mass changed")
        r0=np.hypot(p0["x"].astype(float)-3.,p0["y"].astype(float)-3.)
        bands={}
        for title,lo,hi in (("inner",.5,.75),("disk",.5,2.)):
            take=(r0>=lo)&(r0<hi)
            bands[title]={"radius_range":[lo,hi],
                          "weighted":measure(p0,pw,tnew,take)}
            if not args.weighted_only:
                bands[title]["W0"]=measure(p0,pv,told,take)
                bands[title]["W0_over_weighted_density_abs"]=(
                    bands[title]["W0"]["density_abs_mass_mean"]
                    /bands[title]["weighted"]["density_abs_mass_mean"])
        result["cases"][label]={"time":tnew,
                                "initial_sha256":digest(NEW/label/"initial.dat"),
                                "weighted_endpoint_sha256":digest(new),
                                "weight_range":float(np.ptp(pw["w2"].astype(float))),
                                "bands":bands}
        if not args.weighted_only:
            result["cases"][label]["W0_endpoint_sha256"]=digest(old)
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":main()

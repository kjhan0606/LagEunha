"""Screen a volume-target feedback in the finite-weight Kepler geometry."""

import json
from pathlib import Path
import time

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from export_kepler_q2_weight_rate_csr_t65_20261027 import finite_power_edges
from laguerre_newton_volume_t6_20261026 import newton_jacobian
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from laguerre_weight_rate_actual_t625_20261026 import face_volume_rate
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_weighted_t10_a10_v1")
OUTPUT=Path(__file__).with_name("KEPLER_Q2_WEIGHTED_TARGET_FEEDBACK_T7426_2026-10-27.json")
N=55351
DT=.00193983


def area_metric(area,target,inner):
    error=np.abs(target/area-1.)
    return {"disk_mean_abs_density_equivalent":float(np.mean(error)),
            "disk_p99_abs_density_equivalent":float(np.quantile(error,.99)),
            "inner_mean_abs_density_equivalent":float(np.mean(error[inner])),
            "inner_p99_abs_density_equivalent":float(np.quantile(error[inner],.99))}


def case(label):
    run=ROOT/label
    initial=run/"initial.dat"
    checkpoint=run/"khout.000478.dat"
    t=snapshot_time(checkpoint,N)
    if not 7.42<t<7.43: raise ValueError("Wrong checkpoint")
    p0,_=ordered(initial,N)
    p,_=ordered(checkpoint,N)
    x=np.column_stack((p["x"].astype(float)-3.,p["y"].astype(float)-3.))
    v=np.column_stack((p["vx"].astype(float),p["vy"].astype(float)))
    r0=np.hypot(p0["x"].astype(float)-3.,p0["y"].astype(float)-3.)
    active=np.flatnonzero((r0>=.5)&(r0<2.))
    inner=r0[active]<.75
    lookup=np.full(N,-1,dtype=int)
    lookup[active]=np.arange(len(active))
    W=p["w2"].astype(float)
    area,meta,triangles,centers=power_areas(x,W,active,lookup,return_geometry=True)
    if meta["hidden_active"] or np.any(area<=0): raise ValueError("Bad current diagram")
    target=p0["volume"][active].astype(float)
    fixed,_,nfaces=face_volume_rate(x,v,W,np.zeros(N),active,lookup,triangles,centers)
    J,diagonal=newton_jacobian(x,lookup,len(active),triangles,centers)
    i,j,c1,c2=finite_power_edges(triangles,centers,lookup>=0)
    if len(i)!=nfaces: raise ValueError("Face graph mismatch")
    g=np.linalg.norm(c2-c1,axis=1)/(2*np.linalg.norm(x[i]-x[j],axis=1))
    pair=g*(p["pressure"][i].astype(float)-p["pressure"][j].astype(float))
    h=np.zeros(N)
    np.add.at(h,i,pair)
    np.add.at(h,j,-pair)
    h=h[active]
    h2=float(h@h)
    if h2<=0: raise ValueError("No pressure gradient for projection")
    before=area_metric(area,target,inner)
    candidates={}
    for name,tau in (("zero_volume_rate",None),("feedback_tau_1",1.),
                     ("feedback_tau_0p2",.2)):
        rhs=-fixed if tau is None else -fixed+(target-area)/tau
        rate=np.zeros(len(active))
        for _ in range(16):
            rate+=.8*(rhs-J@rate)/diagonal
        projected=rate-h*float(h@rate)/h2
        allrate=np.zeros(N)
        allrate[active]=projected
        nextW=W+DT*allrate
        nextW-=nextW.min()
        next_area,next_meta=power_areas(x+DT*v,nextW,active,lookup)
        if next_meta["hidden_active"] or np.any(next_area<=0):
            raise ValueError(f"{name} hides a disk power cell")
        guard=float(np.max(np.abs(nextW[i]-nextW[j])/
                           np.sum((x[i]-x[j])**2,axis=1)))
        candidates[name]={
            "tau":tau,"rate_abs_max":float(np.max(np.abs(projected))),
            "projection_relative_change_l2":float(np.linalg.norm(projected-rate)/
                                                   np.linalg.norm(rate)),
            "pressure_weight_power":float(h@projected),
            "jacobi_relative_residual_l2":float(np.linalg.norm(rhs-J@rate)/np.linalg.norm(rhs)),
            "max_weight_difference_over_pair_distance_squared":guard,
            "post_euler":area_metric(next_area,target,inner),
            "fractional_disk_error_change":(
                area_metric(next_area,target,inner)["disk_mean_abs_density_equivalent"]
                /before["disk_mean_abs_density_equivalent"]-1.),
        }
    return {"time":t,"initial_sha256":digest(initial),
            "checkpoint_sha256":digest(checkpoint),
            "active_cells":len(active),"faces_touching_active":nfaces,
            "current":before,"candidates":candidates}


def main():
    started=time.monotonic()
    result={"cases":{label:case(label) for label in ("aligned","jitter")},
            "dt":DT,"elapsed_seconds":time.monotonic()-started,
            "scope":"Read-only one-Euler-step power-geometry screen. A fixed initial-volume target is only a Kepler-equilibrium control; no native RK4 conservation, physical-flow correctness, or long stability is proven."}
    with OUTPUT.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":main()

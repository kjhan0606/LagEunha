"""Measure why native weighted cells still drift by t=7.426, without evolving gas."""

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
OUTPUT=Path(__file__).with_name("KEPLER_Q2_WEIGHTED_JACOBI_RESIDUAL_T7426_2026-10-27.json")
N=55351
DT=.00193983


def values(area,target,rate,inner):
    delta=np.abs(area/target-1.)
    predicted=np.abs(DT*rate/area)
    return {"target_area_mean_abs":float(np.mean(delta)),
            "target_area_p99":float(np.quantile(delta,.99)),
            "inner_target_area_mean_abs":float(np.mean(delta[inner])),
            "predicted_step_area_p99":float(np.quantile(predicted,.99)),
            "predicted_inner_step_area_p99":float(np.quantile(predicted[inner],.99)),
            "rate_l2":float(np.linalg.norm(rate))}


def case(label):
    run=ROOT/label
    initial=run/"initial.dat"
    checkpoint=run/"khout.000478.dat"
    t=snapshot_time(checkpoint,N)
    if not 7.42<t<7.43:
        raise ValueError("Wrong complete native checkpoint time")
    p0,_=ordered(initial,N)
    p,_=ordered(checkpoint,N)
    xy=np.column_stack((p["x"].astype(float)-3.,p["y"].astype(float)-3.))
    xy0=np.column_stack((p0["x"].astype(float)-3.,p0["y"].astype(float)-3.))
    r0=np.linalg.norm(xy0,axis=1)
    active=np.flatnonzero((r0>=.5)&(r0<2.))
    inner=(r0[active]<.75)
    lookup=np.full(N,-1,dtype=int)
    lookup[active]=np.arange(len(active))
    weights=p["w2"].astype(float)
    area,meta,triangles,centers=power_areas(xy,weights,active,lookup,return_geometry=True)
    if meta["hidden_active"] or np.any(area<=0):
        raise ValueError("Hidden weighted power cell")
    native_area=p["volume"][active].astype(float)
    area_error=float(np.max(np.abs(native_area/area-1.)))
    target=p0["volume"][active].astype(float)
    velocity=np.column_stack((p["vx"].astype(float),p["vy"].astype(float)))
    fixed,_,faces=face_volume_rate(xy,velocity,weights,np.zeros(N),active,
                                   lookup,triangles,centers)
    jac,diag=newton_jacobian(xy,lookup,len(active),triangles,centers)
    if np.any(diag<=0):
        raise ValueError("Nonpositive Jacobian diagonal")
    i,j,c1,c2=finite_power_edges(triangles,centers,lookup>=0)
    if len(i)!=faces:
        raise ValueError("Power face graph mismatch")
    g=np.linalg.norm(c2-c1,axis=1)/(2*np.linalg.norm(xy[i]-xy[j],axis=1))
    pressure=p["pressure"].astype(float)
    h=np.zeros(N)
    pair=g*(pressure[i]-pressure[j])
    np.add.at(h,i,pair)
    np.add.at(h,j,-pair)
    h=h[active]
    denominator=float(h@h)
    if denominator<=0: raise ValueError("Zero projected pressure gradient")
    rhs=-fixed
    rate=np.zeros(len(active))
    baseline=values(area,target,fixed,inner)
    history={}
    for iteration in range(1,65):
        rate+=.8*(rhs-jac@rate)/diag
        if iteration in (16,32,64):
            projected=rate-h*(float(h@rate)/denominator)
            residual=fixed+jac@projected
            metrics=values(area,target,residual,inner)
            metrics["jacobi_relative_residual_l2"]=float(
                np.linalg.norm(rhs-jac@rate)/np.linalg.norm(rhs))
            metrics["projection_rate_relative_change_l2"]=float(
                np.linalg.norm(projected-rate)/np.linalg.norm(rate))
            metrics["rate_abs_max"]=float(np.max(np.abs(projected)))
            history[str(iteration)]=metrics
    return {"time":t,"initial_sha256":digest(initial),
            "checkpoint_sha256":digest(checkpoint),"active_cells":len(active),
            "faces_touching_active":faces,"independent_area_error_max":area_error,
            "weight_range":float(np.ptp(weights)),"fixed_weight_rate":baseline,
            "projected_jacobi":history}


def main():
    started=time.monotonic()
    result={"cases":{label:case(label) for label in ("aligned","jitter")},
            "dt_for_predicted_step":DT,"elapsed_seconds":time.monotonic()-started,
            "scope":"Read-only finite-weight geometry/linear residual at native t=7.426. Predicted one-step linear area change is not a native RK4 or long-time stability result."}
    with OUTPUT.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__": main()

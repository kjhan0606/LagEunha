"""Frozen-state pressure-force counterfactual using solved static weights."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_pressure_force_t10_20261027 import (
    finite_faces_touching_active, summarize,
)
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N=55351
GAMMA=float(np.float32(1.666667))


def pressure_force(xy, mass, pressure, active_flag, triangles, centers):
    i,j,c,length=finite_faces_touching_active(triangles,centers,active_flag)
    dr=xy[j]-xy[i]
    distance=np.linalg.norm(dr,axis=1)
    normal=dr/distance[:,None]
    offset=c-.5*(xy[i]+xy[j])
    bi=length[:,None]*(.5*normal+offset/distance[:,None])
    bj=length[:,None]*(.5*normal-offset/distance[:,None])
    dp=pressure[i]-pressure[j]
    force=np.zeros((N,2))
    ia,ja=active_flag[i],active_flag[j]
    np.add.at(force,i[ia],dp[ia,None]*bi[ia])
    np.add.at(force,j[ja],dp[ja,None]*bj[ja])
    acceleration=force/mass[:,None]
    r=np.linalg.norm(xy,axis=1)
    radial=np.einsum("ij,ij->i",xy,acceleration)/np.maximum(r,1e-30)
    tangential=(-xy[:,1]*acceleration[:,0]+xy[:,0]*acceleration[:,1])/np.maximum(r,1e-30)
    return acceleration,radial,tangential,len(i)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("native_receipt",type=Path)
    ap.add_argument("static_weights",type=Path)
    ap.add_argument("native_force_reference",type=Path)
    ap.add_argument("--expected-failed-check",action="append",default=[])
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    receipt=json.loads(args.native_receipt.read_text())
    failed={k for k,v in receipt["checks"].items() if not v}
    if receipt["time"]!=10 or failed!=set(args.expected_failed_check):
        raise ValueError("Unexpected t10 native failed-check label")
    if (digest(Path(receipt["endpoint"]))!=receipt["endpoint_sha256"]
            or digest(Path(receipt["initial"]))!=receipt["initial_sha256"]):
        raise ValueError("Changed native snapshot")
    reference=json.loads(args.native_force_reference.read_text())
    if reference.get("source_snapshot_sha256",reference.get("endpoint_sha256"))!=receipt["endpoint_sha256"]:
        raise ValueError("Native force reference belongs to another state")
    p0,_=ordered(Path(receipt["initial"]),N)
    p,_=ordered(Path(receipt["endpoint"]),N)
    xy0=np.column_stack((p0["x"].astype(float)-3.,p0["y"].astype(float)-3.))
    xy=np.column_stack((p["x"].astype(float)-3.,p["y"].astype(float)-3.))
    r0=np.linalg.norm(xy0,axis=1)
    active=np.flatnonzero((r0>=.5)&(r0<2.))
    active_flag=np.zeros(N,dtype=bool)
    active_flag[active]=True
    lookup=np.full(N,-1,dtype=int)
    lookup[active]=np.arange(len(active))
    with np.load(args.static_weights) as saved:
        if not np.array_equal(saved["active"],active):
            raise ValueError("Wrong static-weight active IDs")
        local_weight=saved["weight"].copy()
    weight=np.zeros(N)
    weight[active]=local_weight
    native_area=p["volume"].astype(float)
    mass=p["mass"].astype(float)
    native_pressure=p["pressure"].astype(float)
    w0_area,w0_meta,w0_triangles,w0_centers=power_areas(
        xy,np.zeros(N),active,lookup,return_geometry=True)
    if w0_meta["hidden_active"] or np.max(np.abs(w0_area/native_area[active]-1))>1e-5:
        raise ValueError("W0 geometry differs from native area")
    f0,r0f,t0f,_=pressure_force(xy,mass,native_pressure,active_flag,
                                 w0_triangles,w0_centers)
    w0_summary=summarize(active_flag,mass,native_pressure,f0[:,0],f0[:,1],r0f,t0f)
    if abs(w0_summary["acceleration_mass_weighted_rms"]-
           reference["whole_disk"]["acceleration_mass_weighted_rms"])>1e-12:
        raise ValueError("Independent W0 force reconstruction disagrees with reference")
    area,meta,triangles,centers=power_areas(
        xy,weight,active,lookup,return_geometry=True)
    if (meta["hidden_active"] or np.any(area<=0)
            or np.max(np.abs(area/p0["volume"][active]-1))>1e-7):
        raise ValueError("Static target-volume power diagram invalid")
    i,j,_,_=finite_faces_touching_active(triangles,centers,active_flag)
    adjacent=np.unique(np.r_[i[~active_flag[i]],j[~active_flag[j]]])
    selected=np.unique(np.r_[active,adjacent])
    union_lookup=np.full(N,-1,dtype=int)
    union_lookup[selected]=np.arange(len(selected))
    selected_area,selected_meta=power_areas(xy,weight,selected,union_lookup)
    if selected_meta["hidden_active"] or np.any(selected_area<=0):
        raise ValueError("Hidden or invalid disk-adjacent weighted cell")
    if np.max(np.abs(selected_area[union_lookup[active]]/area-1))>1e-9:
        raise ValueError("Power area changed while adding inactive neighbors")
    adiabatic_pressure=native_pressure.copy()
    adiabatic_pressure[selected]*=(native_area[selected]/selected_area)**GAMMA
    f_static,r_static,t_static,_=pressure_force(
        xy,mass,native_pressure,active_flag,triangles,centers)
    f_adiabatic,r_adiabatic,t_adiabatic,faces=pressure_force(
        xy,mass,adiabatic_pressure,active_flag,triangles,centers)
    mid=active_flag&(r0>=.75)&(r0<1.75)
    def pack(acc,rad,tan,pressure):
        return {"whole_disk":summarize(active_flag,mass,pressure,acc[:,0],acc[:,1],rad,tan),
                "middle_disk_0p75_to_1p75":summarize(mid,mass,pressure,acc[:,0],acc[:,1],rad,tan)}
    result={"time":10.0,"native_receipt":str(args.native_receipt.resolve()),
            "endpoint_sha256":receipt["endpoint_sha256"],
            "native_failed_checks":sorted(failed),
            "weights_sha256":digest(args.static_weights),
            "native_force_reference_sha256":digest(args.native_force_reference),
            "active_disk_sites":len(active),"adjacent_inactive_sites":len(adjacent),
            "weighted_faces_touching_disk":faces,
            "target_area_max_relative_error":float(np.max(np.abs(area/p0["volume"][active]-1))),
            "pressure_reset":"Adiabatic P_new=P_native*(V_native/V_power)^gamma for every disk and directly adjacent inactive cell; mass and entropy labels held fixed.",
            "w0_native_pressure":pack(f0,r0f,t0f,native_pressure),
            "static_weight_native_pressure":pack(f_static,r_static,t_static,native_pressure),
            "static_weight_adiabatic_pressure":pack(f_adiabatic,r_adiabatic,t_adiabatic,adiabatic_pressure),
            "scope":"Frozen-state counterfactual force at fixed sites. Reweighted EOS fields have not been evolved, coupled to RK4 or checked for generalized weight-work conservation. Native failed gate remains failed."}
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    for key in ("w0_native_pressure","static_weight_native_pressure",
                "static_weight_adiabatic_pressure"):
        print(key,result[key]["whole_disk"]["acceleration_mass_weighted_rms"])


if __name__=="__main__":
    main()

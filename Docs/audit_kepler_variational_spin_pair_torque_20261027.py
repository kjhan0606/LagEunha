"""Measure a reciprocal face-force decomposition of native variational pressure.

The zero-sum face gauge converts the directed variational force into an
antisymmetric pair force without changing any closed-cell net force. Its
moment about each generator is a spin-torque diagnostic, not a unique
physical traction or an implementation of cell spin.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N=55351
BANDS=((.5,.75),(.75,1.),(1.,1.25),(1.25,1.5),(1.5,1.75),(1.75,2.))
EPS=.01175


def cross(a,b):
    return a[:,0]*b[:,1]-a[:,1]*b[:,0]


def stats(mask,mass,torque,face_torque_abs,scale,dt):
    values=np.abs(torque[mask])
    fractions=dt*values/scale[mask]
    m=mass[mask]
    face_abs=float(np.sum(face_torque_abs[mask],dtype=np.longdouble))
    net_abs=float(np.sum(values,dtype=np.longdouble))
    return {"cells":int(mask.sum()),
            "signed_torque_sum":float(np.sum(torque[mask],dtype=np.longdouble)),
            "absolute_torque_sum":net_abs,
            "absolute_face_torque_sum":face_abs,
            "net_cell_absolute_over_face_absolute":net_abs/max(face_abs,1e-300),
            "abs_torque_per_mass_rms":float(np.sqrt(np.average((values/m)**2,weights=m))),
            "dt_abs_torque_over_m_area_omega_mass_mean":float(np.average(fractions,weights=m)),
            "dt_abs_torque_over_m_area_omega_p95":float(np.quantile(fractions,.95)),
            "dt_abs_torque_over_m_area_omega_p99":float(np.quantile(fractions,.99)),
            "dt_abs_torque_over_m_area_omega_max":float(np.max(fractions))}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("receipt",type=Path)
    ap.add_argument("--expected-failed-check",action="append",default=[])
    ap.add_argument("--gpu-reference",type=Path)
    ap.add_argument("--gpu-case",choices=("aligned","jitter"))
    ap.add_argument("--dt",type=float,default=.00193983)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    if not np.isfinite(args.dt) or args.dt<=0:
        raise ValueError("Expected positive finite diagnostic dt")
    receipt=json.loads(args.receipt.read_text())
    failed={k for k,v in receipt["checks"].items() if not v}
    if failed!=set(args.expected_failed_check):
        raise ValueError(f"Unexpected receipt failures: {sorted(failed)}")
    for key in ("initial","endpoint"):
        if digest(Path(receipt[key]))!=receipt[key+"_sha256"]:
            raise ValueError(f"Changed {key} snapshot")
    initial,_=ordered(Path(receipt["initial"]),N)
    p,meta_snapshot=ordered(Path(receipt["endpoint"]),N)
    if abs(float(meta_snapshot["t"])-float(receipt["time"]))>1e-10:
        raise ValueError("Endpoint time differs from receipt")
    if not np.array_equal(initial["mass"],p["mass"]):
        raise ValueError("Expected this experiment's exact fixed cell masses")
    xy0=np.column_stack((initial["x"].astype(float)-3.,initial["y"].astype(float)-3.))
    xy=np.column_stack((p["x"].astype(float)-3.,p["y"].astype(float)-3.))
    r0=np.linalg.norm(xy0,axis=1)
    active=np.flatnonzero((r0>=.5)&(r0<2.))
    active_flag=np.zeros(N,dtype=bool)
    active_flag[active]=True
    lookup=np.full(N,-1,dtype=int)
    lookup[active]=np.arange(len(active))
    weight=p["w2"].astype(float)
    area,meta,triangles,centers=power_areas(xy,weight,active,lookup,return_geometry=True)
    if meta["hidden_active"] or np.any(area<=0):
        raise ValueError("Invalid active power-cell geometry")
    area_error=float(np.max(np.abs(area/p["volume"][active].astype(float)-1.)))
    if area_error>1e-5:
        raise ValueError("Independent power areas differ from native snapshot")
    i,j,c,length=finite_faces_touching_active(triangles,centers,active_flag)
    dr=xy[j]-xy[i]
    d=np.linalg.norm(dr,axis=1)
    if np.any(d<=0) or np.any(length<=0):
        raise ValueError("Degenerate face or site pair")
    n=dr/d[:,None]
    area_vector=length[:,None]*n
    offset=c-.5*(xy[i]+xy[j])
    bi=length[:,None]*(.5*n+offset/d[:,None])
    bj=length[:,None]*(.5*n-offset/d[:,None])
    pressure=p["pressure"].astype(float)
    mass=p["mass"].astype(float)
    dp=pressure[i]-pressure[j]
    directed_i=dp[:,None]*bi
    directed_j=dp[:,None]*bj
    pair_force=directed_i-pressure[i,None]*area_vector
    pair_force_j=directed_j+pressure[j,None]*area_vector
    pair_residual=pair_force+pair_force_j
    native_force=np.zeros((N,2))
    pair_force_sum=np.zeros((N,2))
    closure=np.zeros((N,2))
    spin_torque=np.zeros(N)
    face_torque_abs=np.zeros(N)
    moment_closure=np.zeros(N)
    ia,ja=active_flag[i],active_flag[j]
    np.add.at(native_force,i[ia],directed_i[ia])
    np.add.at(native_force,j[ja],directed_j[ja])
    np.add.at(pair_force_sum,i[ia],pair_force[ia])
    np.add.at(pair_force_sum,j[ja],-pair_force[ja])
    np.add.at(closure,i[ia],area_vector[ia])
    np.add.at(closure,j[ja],-area_vector[ja])
    ti=cross(c-xy[i],pair_force)
    tj=cross(c-xy[j],-pair_force)
    orbital=cross(xy[i]-xy[j],pair_force)
    np.add.at(spin_torque,i[ia],ti[ia])
    np.add.at(spin_torque,j[ja],tj[ja])
    np.add.at(face_torque_abs,i[ia],np.abs(ti[ia]))
    np.add.at(face_torque_abs,j[ja],np.abs(tj[ja]))
    moment_i=cross(c-xy[i],area_vector)
    moment_j=cross(c-xy[j],-area_vector)
    np.add.at(moment_closure,i[ia],moment_i[ia])
    np.add.at(moment_closure,j[ja],moment_j[ja])
    face_torque_formula_error=max(float(np.max(np.abs(ti+pressure[i]*moment_i))),
                                  float(np.max(np.abs(tj+pressure[j]*moment_j))))
    torque_residual=orbital+ti+tj
    pair_scale=max(float(np.max(np.linalg.norm(pair_force,axis=1))),1e-300)
    torque_scale=max(float(np.max(np.abs(ti)+np.abs(tj)+np.abs(orbital))),1e-300)
    gauge_accel_error=float(np.max(np.linalg.norm(
        (native_force[active]-pair_force_sum[active])/mass[active,None],axis=1)))
    closure_max=float(np.max(np.linalg.norm(closure[active],axis=1)))
    pair_force_rel=float(np.max(np.linalg.norm(pair_residual,axis=1))/pair_scale)
    pair_torque_rel=float(np.max(np.abs(torque_residual))/torque_scale)
    face_torque_scale=max(float(np.max(np.abs(ti)+np.abs(tj))),1e-300)
    face_formula_rel=face_torque_formula_error/face_torque_scale
    net_face_ratio=float(np.sum(np.abs(spin_torque[active]),dtype=np.longdouble)
                         /max(float(np.sum(face_torque_abs[active],dtype=np.longdouble)),1e-300))
    moment_closure_max=float(np.max(np.abs(moment_closure[active])))
    omega=(np.sum(xy*xy,axis=1)+EPS*EPS)**(-.75)
    scale=mass*p["volume"].astype(float)*omega
    if np.any(scale[active]<=0) or not np.all(np.isfinite(scale[active])):
        raise ValueError("Invalid cell angular-momentum scale")
    bands=[]
    for lo,hi in BANDS:
        mask=active_flag&(r0>=lo)&(r0<hi)
        bands.append({"r_min":lo,"r_max":hi,
                      **stats(mask,mass,spin_torque,face_torque_abs,scale,args.dt)})
    gpu_check=None
    if args.gpu_reference:
        if not args.gpu_case or abs(float(receipt["time"])-6.5)>1e-10 or np.any(weight!=0):
            raise ValueError("GPU sample reference requires explicit W0 t=6.5 case")
        capture=json.loads(args.gpu_reference.read_text())["captures"][args.gpu_case]["cells"]
        errors=[]
        for key,row in capture.items():
            idx=int(key)
            if not active_flag[idx]:
                raise ValueError("GPU sample is outside disk")
            expected=np.array([row["pressure_acceleration_x"],row["pressure_acceleration_y"]])
            errors.append(float(np.linalg.norm(native_force[idx]/mass[idx]-expected)))
        gpu_check={"case":args.gpu_case,"cells":len(errors),
                   "reference_sha256":digest(args.gpu_reference),
                   "max_acceleration_vector_difference":max(errors)}
    checks={
        "independent_area":area_error<1e-5,
        "closed_disk_cells":closure_max<1e-9,
        "reciprocal_pair_force":pair_force_rel<1e-12,
        "same_native_net_force":gauge_accel_error<1e-10,
        "orbital_plus_spin_pair_torque":pair_torque_rel<1e-12,
        "face_torque_identity":face_formula_rel<1e-12,
        "closed_cell_moment_identity":moment_closure_max<1e-10 and net_face_ratio<1e-10,
        "native_gpu_sample":gpu_check is None or gpu_check["max_acceleration_vector_difference"]<5e-11,
    }
    result={"receipt":str(args.receipt.resolve()),"receipt_sha256":digest(args.receipt),
            "initial_sha256":receipt["initial_sha256"],"endpoint_sha256":receipt["endpoint_sha256"],
            "time":float(receipt["time"]),"native_failed_checks":sorted(failed),
            "diagnostic_dt":args.dt,"disk_sites":len(active),"faces_touching_disk":len(i),
            "disk_disk_faces":int(np.count_nonzero(ia&ja)),"weighted":bool(np.any(weight!=0)),
            "independent_power_area_error_max":area_error,"cell_area_vector_closure_max":closure_max,
            "reciprocal_pair_force_relative_max":pair_force_rel,
            "gauge_vs_native_acceleration_max":gauge_accel_error,
            "orbital_plus_spin_pair_torque_relative_max":pair_torque_rel,
            "face_torque_identity_relative_max":face_formula_rel,
            "cell_moment_closure_max":moment_closure_max,
            "gpu_sample_comparison":gpu_check,
            "whole_disk":stats(active_flag,mass,spin_torque,face_torque_abs,scale,args.dt),
            "radial_bands":bands,"checks":checks,
            "scope":"Frozen-state reconstructed variational pressure. The chosen reciprocal "
                    "face-force gauge preserves closed-cell net force and pair angular momentum, "
                    "but face torque allocation is not unique; no native spin state, energy "
                    "exchange, Riemann flux or time-integrated spin is represented."}
    with args.output.open("x") as output:
        json.dump(result,output,indent=2,allow_nan=False)
        output.write("\n")
    print(json.dumps({"time":result["time"],"weighted":result["weighted"],
                      "failed_checks":[k for k,v in checks.items() if not v],
                      "disk":result["whole_disk"]},indent=2))
    if not all(checks.values()):
        raise SystemExit("Variational pair-torque reconstruction failed")


if __name__=="__main__":
    main()

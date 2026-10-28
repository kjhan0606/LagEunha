"""Reconstruct the native W=0 variational pressure force from saved faces."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
BANDS = ((.5,.75),(.75,1.),(1.,1.25),(1.25,1.5),(1.5,1.75),(1.75,2.))


def finite_faces_touching_active(triangles, centers, active_flag):
    edges = np.vstack((triangles[:,[0,1]],triangles[:,[1,2]],
                       triangles[:,[2,0]]))
    tids = np.tile(np.arange(len(triangles)),3)
    edges.sort(axis=1)
    order = np.lexsort((edges[:,1],edges[:,0]))
    edges,tids = edges[order],tids[order]
    starts = np.r_[0,1+np.flatnonzero(np.any(edges[1:] != edges[:-1],axis=1))]
    count = np.diff(np.r_[starts,len(edges)])
    touch = active_flag[edges[starts,0]] | active_flag[edges[starts,1]]
    if np.any(count[touch] != 2):
        raise ValueError("Non-finite power face touches disk")
    starts = starts[touch]
    i,j = edges[starts,0],edges[starts,1]
    c1,c2=centers[tids[starts]],centers[tids[starts+1]]
    return i,j,.5*(c1+c2),np.linalg.norm(c2-c1,axis=1)


def weighted_rms(values,mass):
    return float(np.sqrt(np.sum(mass*values*values,dtype=np.longdouble)
                         /np.sum(mass,dtype=np.longdouble)))


def summarize(mask,mass,pressure,ax,ay,ar,at):
    magnitude=np.hypot(ax[mask],ay[mask])
    p=pressure[mask]
    return {"cells":int(np.count_nonzero(mask)),
            "pressure_mean":float(np.mean(p)),
            "pressure_relative_rms":float(np.std(p)/np.mean(p)),
            "acceleration_median":float(np.median(magnitude)),
            "acceleration_p95":float(np.quantile(magnitude,.95)),
            "acceleration_p99":float(np.quantile(magnitude,.99)),
            "acceleration_mass_weighted_rms":weighted_rms(magnitude,mass[mask]),
            "radial_acceleration_mass_weighted_rms":weighted_rms(ar[mask],mass[mask]),
            "tangential_acceleration_mass_weighted_rms":weighted_rms(at[mask],mass[mask])}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("native_receipt",type=Path)
    ap.add_argument("--checkpoint",type=Path,
                    help="Complete intermediate khout in the same native run; default is receipt endpoint")
    ap.add_argument("--initial-state",action="store_true",
                    help="Read the hash-verified native receipt initial snapshot at t=0")
    ap.add_argument("--expected-failed-check",action="append",default=[])
    ap.add_argument("--gpu_reference",type=Path)
    ap.add_argument("--gpu_case",choices=("aligned","jitter"))
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    if args.initial_state and args.checkpoint:
        ap.error("--initial-state and --checkpoint are mutually exclusive")
    receipt=json.loads(args.native_receipt.read_text())
    failed={k for k,v in receipt["checks"].items() if not v}
    if (failed != set(args.expected_failed_check) or
            receipt["time"] not in (6.5,10.0,20.0,30.0)):
        raise ValueError(f"Unexpected failed checks or time: {sorted(failed)}")
    for key,hash_key in (("initial","initial_sha256"),
                         ("endpoint","endpoint_sha256")):
        if digest(Path(receipt[key])) != receipt[hash_key]:
            raise ValueError(f"Changed native {key} snapshot")
    p0,_=ordered(Path(receipt["initial"]),N)
    checkpoint=(Path(receipt["initial"]) if args.initial_state else
                args.checkpoint or Path(receipt["endpoint"]))
    if not args.initial_state and checkpoint.resolve().parent != Path(receipt["endpoint"]).resolve().parent:
        raise ValueError("Intermediate checkpoint is not in the audited native run")
    p,checkpoint_meta=ordered(checkpoint,N)
    checkpoint_time=float(checkpoint_meta["t"])
    if args.initial_state and abs(checkpoint_time)>1e-12:
        raise ValueError("Expected initial snapshot at t=0")
    if not args.initial_state and not (0.0 <= checkpoint_time <= float(receipt["time"])+1e-10):
        raise ValueError("Intermediate checkpoint is outside the audited segment")
    checkpoint_sha=digest(checkpoint)
    if args.initial_state and checkpoint_sha != receipt["initial_sha256"]:
        raise ValueError("Initial SHA differs from native receipt")
    if not args.initial_state and not args.checkpoint and checkpoint_sha != receipt["endpoint_sha256"]:
        raise ValueError("Default endpoint SHA differs from native receipt")
    if not np.array_equal(p0["mass"],p["mass"]) or np.any(p["w2"] != 0):
        raise ValueError("Expected fixed mass and W=0")
    xy0=np.column_stack((p0["x"].astype(float)-3.,p0["y"].astype(float)-3.))
    xy=np.column_stack((p["x"].astype(float)-3.,p["y"].astype(float)-3.))
    r0=np.linalg.norm(xy0,axis=1)
    active=np.flatnonzero((r0>=.5)&(r0<2.))
    active_flag=np.zeros(N,dtype=bool)
    active_flag[active]=True
    lookup=np.full(N,-1,dtype=int)
    lookup[active]=np.arange(len(active))
    area,meta,triangles,centers=power_areas(
        xy,np.zeros(N),active,lookup,return_geometry=True)
    if meta["hidden_active"] or np.any(area<=0):
        raise ValueError("Invalid W=0 reference Voronoi geometry")
    area_error=np.max(np.abs(area/p["volume"][active]-1))
    if area_error>1e-5:
        raise ValueError("Reference Voronoi area differs from native saved area")
    i,j,c,length=finite_faces_touching_active(triangles,centers,active_flag)
    dr=xy[j]-xy[i]
    distance=np.linalg.norm(dr,axis=1)
    normal=dr/distance[:,None]
    offset=c-.5*(xy[i]+xy[j])
    bi=length[:,None]*(.5*normal+offset/distance[:,None])
    bj=length[:,None]*(.5*normal-offset/distance[:,None])
    pressure=p["pressure"].astype(float)
    dp=pressure[i]-pressure[j]
    force=np.zeros((N,2))
    ia,ja=active_flag[i],active_flag[j]
    np.add.at(force,i[ia],dp[ia,None]*bi[ia])
    np.add.at(force,j[ja],dp[ja,None]*bj[ja])
    mass=p["mass"].astype(float)
    acceleration=force/mass[:,None]
    rr=np.linalg.norm(xy,axis=1)
    radial=np.einsum("ij,ij->i",xy,acceleration)/np.maximum(rr,1e-30)
    tangent=(-xy[:,1]*acceleration[:,0]+xy[:,0]*acceleration[:,1])/np.maximum(rr,1e-30)
    disk=summarize(active_flag,mass,pressure,acceleration[:,0],acceleration[:,1],radial,tangent)
    bands=[]
    for lo,hi in BANDS:
        mask=active_flag&(r0>=lo)&(r0<hi)
        bands.append({"r_min":lo,"r_max":hi,
                      **summarize(mask,mass,pressure,acceleration[:,0],acceleration[:,1],radial,tangent)})
    result={"native_receipt":str(args.native_receipt.resolve()),
            "native_receipt_sha256":digest(args.native_receipt),
            "source_snapshot":str(checkpoint.resolve()),
            "source_snapshot_sha256":checkpoint_sha,
            "receipt_endpoint_sha256":receipt["endpoint_sha256"],
            "receipt_endpoint_time":float(receipt["time"]),
            "time":checkpoint_time,
            "initial_state":args.initial_state,
            "native_failed_checks":sorted(failed),
            "disk_sites":len(active),"faces_touching_disk":len(i),
            "reference_area_max_relative_difference_from_native":float(area_error),
            "force_formula":"Each i<j finite face: F_i += (P_i-P_j)L[n/2+(c-mid)/d]; F_j += (P_i-P_j)L[n/2-(c-mid)/d]. Native W=0 variational pressure only.",
            "whole_disk":disk,"radial_bands":bands,
            "scope":"Frozen-state pressure-force reconstruction from saved cell pressure and independent Voronoi faces. Native failed-check label belongs to the audited receipt endpoint; an initial/intermediate snapshot is exploratory and has no independent full native gate. No RK4 force history, gravity, or native all-cell pressure dump."}
    if args.gpu_reference:
        if not args.gpu_case or checkpoint_time != 6.5:
            raise ValueError("GPU reference needs t=6.5 and explicit case")
        gpu=json.loads(args.gpu_reference.read_text())["captures"][args.gpu_case]["cells"]
        errs=[]
        for key,row in gpu.items():
            idx=int(key)
            if not active_flag[idx]:
                raise ValueError("GPU sample is outside disk")
            err=np.linalg.norm(acceleration[idx]-[row["pressure_acceleration_x"],
                                                  row["pressure_acceleration_y"]])
            errs.append(err)
        result["gpu_sample_comparison"]={
            "reference_sha256":digest(args.gpu_reference),
            "case":args.gpu_case,"cells":len(errs),
            "median_absolute_acceleration_vector_difference":float(np.median(errs)),
            "p95_absolute_acceleration_vector_difference":float(np.quantile(errs,.95)),
            "max_absolute_acceleration_vector_difference":float(np.max(errs))}
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps({"time":result["time"],"failed":result["native_failed_checks"],
                      "disk_accel_rms":disk["acceleration_mass_weighted_rms"],
                      "gpu_compare":result.get("gpu_sample_comparison")}))


if __name__=="__main__":
    main()

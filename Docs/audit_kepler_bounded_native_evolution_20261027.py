"""Audit a native fixed-polygon Kepler endpoint against exact circular sites."""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_native_ballistic_20261027 import disk_geometry
from audit_kepler_sparse_state_consistency_20261027 import state_metrics
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N = 55351
DOMAIN_AREA = .5*1024*2.75**2*np.sin(2*np.pi/1024)
APOTHEM = 2.75*np.cos(np.pi/1024)
GAMMA_MINUS_ONE = float(np.float32(1.666667))-1.


def energy(p):
    x=p["x"].astype(float)-3.
    y=p["y"].astype(float)-3.
    mass=p["mass"].astype(float)
    vx=p["vx"].astype(float)
    vy=p["vy"].astype(float)
    return float(np.sum(.5*mass*(vx*vx+vy*vy)+p["ie"].astype(float)
                        -mass/np.sqrt(x*x+y*y+.01175**2),dtype=np.longdouble))


def weighted_abs(values,weights):
    return float(np.sum(np.abs(values)*weights,dtype=np.longdouble)
                 /np.sum(weights,dtype=np.longdouble))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("initial",type=Path)
    ap.add_argument("run",type=Path)
    ap.add_argument("--time",type=float,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--through-endpoint",action="store_true",
                    help="Limit log checks to the completed endpoint step marker")
    ap.add_argument("--expected-binary-sha256",
                    help="Also check the native executable hash")
    args=ap.parse_args()
    if abs(snapshot_time(args.initial,N))>1e-12:
        raise ValueError("Reference is not complete exact native t=0")
    candidates=[path for path in args.run.glob("khout.*.dat")
                if abs(snapshot_time(path,N)-args.time)<1e-10]
    if len(candidates)!=1:
        raise ValueError(f"Expected exactly one complete t={args.time} endpoint")
    endpoint=candidates[0]
    p0,_=ordered(args.initial,N)
    p,_=ordered(endpoint,N)
    s0,s1=state_metrics(p0),state_metrics(p)
    x0=p0["x"].astype(float)-3.
    y0=p0["y"].astype(float)-3.
    r0=np.hypot(x0,y0)
    disk=(r0>=.5)&(r0<2.)
    outer=(r0>=2.35)
    selected=np.flatnonzero(disk)
    xy,ballistic_area,image_safety=disk_geometry(p0,args.time,selected)
    initial_density=p0["den"][selected].astype(float)
    density=p["den"][selected].astype(float)
    mass=p0["mass"][selected].astype(float)
    volume=p["volume"][selected].astype(float)
    native_change=density/initial_density-1.
    ballistic_change=p0["volume"][selected].astype(float)/ballistic_area-1.
    pointwise_residual=ballistic_area/volume-1.
    position_error=np.hypot(p["x"][selected].astype(float)-xy[selected,0],
                            p["y"][selected].astype(float)-xy[selected,1])
    entropy_log=(np.log(p["ie"][selected].astype(float)/p0["ie"][selected].astype(float))
                 +GAMMA_MINUS_ONE*np.log(volume/p0["volume"][selected].astype(float)))
    x=p["x"].astype(float)-3.
    y=p["y"].astype(float)-3.
    r=np.hypot(x,y)
    vx=p["vx"].astype(float)
    vy=p["vy"].astype(float)
    vr=(x*vx+y*vy)/np.maximum(r,1e-30)
    vk=r0*(r0*r0+.01175**2)**(-.75)
    disk_vr_rms=float(np.sqrt(np.sum(mass*(vr[disk]/vk[disk])**2,dtype=np.longdouble)
                              /np.sum(mass,dtype=np.longdouble)))
    L0=float(np.sum(mass*(x0[disk]*p0["vy"][disk].astype(float)
                          -y0[disk]*p0["vx"][disk].astype(float)),dtype=np.longdouble))
    L1=float(np.sum(mass*(x[disk]*vy[disk]-y[disk]*vx[disk]),dtype=np.longdouble))
    e0,e1=energy(p0),energy(p)
    log=(args.run/"log").read_text(errors="replace")
    if args.through_endpoint:
        step=int(endpoint.name.split(".")[1])
        log_lines=log.splitlines()
        marker=next((i for i,line in enumerate(log_lines)
                     if line.startswith("Time=") and
                     re.search(rf"\bstep=\s*{step}\b",line)),None)
        if marker is None:
            raise ValueError("Missing complete log marker for endpoint step")
        log="\n".join(log_lines[:marker+1])
    lines=[line for line in log.splitlines() if "[RK4E]" in line]
    if not lines:
        raise ValueError("Missing native RK4 energy ledger")
    ledger={}
    for key in ("floor_cum","de_cum","sfl_cum"):
        match=re.search(rf"\b{key}=([0-9.eE+-]+)",lines[-1])
        if not match:
            raise ValueError(f"Missing {key}")
        ledger[key]=float(match.group(1))
    density_abs=weighted_abs(native_change,mass)
    ballistic_abs=weighted_abs(ballistic_change,mass)
    residual_abs=weighted_abs(pointwise_residual,mass)
    entropy_abs=weighted_abs(entropy_log,mass)
    position_abs=weighted_abs(position_error,mass)
    energy_rel=abs(e1-e0)/abs(e0)
    angular_drift=L1/L0-1.
    checks={
        "finite_positive":bool(s0["positive"] and s1["positive"]),
        "fixed_mass_exact":bool(np.array_equal(p0["mass"],p["mass"])),
        "W0_exact":bool(np.all(p0["w2"]==0) and np.all(p["w2"]==0)),
        "initial_area_tiles_polygon":abs(s0["volume_sum"]-DOMAIN_AREA)<1e-6,
        "endpoint_area_tiles_polygon":abs(s1["volume_sum"]-DOMAIN_AREA)<1e-5,
        "density_and_EOS":max(s0["density_relative_max"],s0["eos_pressure_relative_max"],
                              s1["density_relative_max"],s1["eos_pressure_relative_max"])<1e-6,
        "zero_injection":all(value==0.0 for value in ledger.values()),
        "total_energy":energy_rel<1e-3,
        "disk_density":density_abs<.05,
        "disk_entropy":entropy_abs<.01,
        "disk_orbit":position_abs<.01,
        "disk_ballistic_residual":residual_abs<1e-3,
        "disk_radial_velocity":disk_vr_rms<.05,
        "disk_angular_momentum":abs(angular_drift)<.01,
        "inside_polygon_incircle":bool(np.all(r<APOTHEM)),
        "no_wall_crossing":"[FIXED_POLYGON_CROSSING]" not in log,
        "fixed_polygon_enabled":"[FIXED_POLYGON] enabled" in log,
    }
    if args.through_endpoint:
        checks["no_geometry_abort"]=("[GPU_TESS_INVALID]" not in log and
                                     "[GPU_CLIP_CAPACITY]" not in log)
    if args.expected_binary_sha256:
        checks["expected_binary"]=(digest(args.run/"eunha2.exe") ==
                                   args.expected_binary_sha256)
    checks={key:bool(value) for key,value in checks.items()}
    result={
        "initial":str(args.initial.resolve()),"endpoint":str(endpoint.resolve()),
        "initial_sha256":digest(args.initial),"endpoint_sha256":digest(endpoint),
        "time":args.time,"n":N,"disk_sites":int(np.sum(disk)),
        "initial_state":s0,"endpoint_state":s1,
        "total_energy_relative":energy_rel,"rk4_ledger":ledger,
        "disk_density_abs_mass_weighted":density_abs,
        "disk_ballistic_density_abs_mass_weighted":ballistic_abs,
        "disk_native_minus_ballistic_abs_mass_weighted":residual_abs,
        "disk_entropy_log_abs_mass_weighted":entropy_abs,
        "disk_position_error_abs_mass_weighted":position_abs,
        "disk_radial_speed_rms_over_circular":disk_vr_rms,
        "disk_angular_momentum_relative":angular_drift,
        "outer_density_abs_mass_weighted":weighted_abs(
            p["den"][outer].astype(float)/p0["den"][outer].astype(float)-1.,
            p0["mass"][outer].astype(float)),
        "minimum_ballistic_image_safety_ratio":image_safety,
        "max_current_site_radius":float(np.max(r)),
        "checks":checks,
        "scope":"One-rank native fixed-polygon CUDA W0 endpoint; requires later boundary, scaling, and long-time validation",
    }
    if args.through_endpoint:
        result["log_scope"]="Completed steps through endpoint marker only"
    if args.expected_binary_sha256:
        result["binary_sha256"]=digest(args.run/"eunha2.exe")
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2)
        stream.write("\n")
    print(json.dumps(result,indent=2))
    if not all(checks.values()):
        raise SystemExit("Bounded native Kepler evolution gate failed")


if __name__=="__main__":
    main()

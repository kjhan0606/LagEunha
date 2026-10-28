"""Gate a single native weighted Kepler RK4 replay from accepted t=6.5."""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_state_consistency_20261027 import state_metrics
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N = 55351
AREA = .5*1024*2.75**2*np.sin(2*np.pi/1024)


def energy(p):
    x=p["x"].astype(float)-3.
    y=p["y"].astype(float)-3.
    m=p["mass"].astype(float)
    vx=p["vx"].astype(float)
    vy=p["vy"].astype(float)
    return float(np.sum(.5*m*(vx*vx+vy*vy)+p["ie"].astype(float)
                        -m/np.sqrt(x*x+y*y+.01175**2), dtype=np.longdouble))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("run",type=Path)
    parser.add_argument("--time",type=float,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    restart=args.run/"khout.000001.dat"
    t0=snapshot_time(restart,N)
    if abs(t0-6.5)>1e-10 or not args.time>t0:
        raise ValueError("Expected accepted t=6.5 restart")
    endpoints=[path for path in args.run.glob("khout.*.dat")
               if abs(snapshot_time(path,N)-args.time)<1e-10]
    if len(endpoints)!=1:
        raise ValueError("Expected one complete weighted endpoint")
    endpoint=endpoints[0]
    p0,_=ordered(restart,N)
    p1,_=ordered(endpoint,N)
    s0,s1=state_metrics(p0),state_metrics(p1)
    xy0=np.column_stack((p0["x"].astype(float)-3.,p0["y"].astype(float)-3.))
    xy1=np.column_stack((p1["x"].astype(float)-3.,p1["y"].astype(float)-3.))
    r0=np.linalg.norm(xy0,axis=1)
    active=(r0>=.5)&(r0<2.)
    ids=np.flatnonzero(active)
    lookup=np.full(N,-1,dtype=int)
    lookup[ids]=np.arange(len(ids))
    weight=p1["w2"].astype(float)
    if not np.all(np.isfinite(weight)):
        raise ValueError("Nonfinite weighted endpoint")
    independent_area,meta=power_areas(xy1,weight,ids,lookup)
    if meta["hidden_active"] or np.any(independent_area<=0):
        raise ValueError("Hidden weighted disk power cell")
    native_area=p1["volume"][ids].astype(float)
    area_error=np.abs(native_area/independent_area-1.)
    phase=(args.time-t0)*(r0*r0+.01175**2)**(-.75)
    x_exact=xy0[:,0]*np.cos(phase)-xy0[:,1]*np.sin(phase)
    y_exact=xy0[:,0]*np.sin(phase)+xy0[:,1]*np.cos(phase)
    orbit_max=float(np.max(np.hypot(xy1[ids,0]-x_exact[ids],
                                      xy1[ids,1]-y_exact[ids])))
    e0,e1=energy(p0),energy(p1)
    energy_relative=abs(e1-e0)/abs(e0)
    log=(args.run/"log").read_text(errors="replace")
    stages=[]
    pattern=re.compile(r"\[GPU_VOLW2\] stage=(\d+) active=(\d+) faces=(\d+) "
                       r"alpha=([\deE+.-]+) work=([\deE+.-]+) elapsed=([\deE+.-]+) s")
    for line in log.splitlines():
        match=pattern.fullmatch(line)
        if match:
            stage,nactive,faces,alpha,work,seconds=match.groups()
            stages.append({"stage":int(stage),"active":int(nactive),
                           "faces":int(faces),"alpha":float(alpha),
                           "projected_pressure_weight_power":float(work),
                           "elapsed_seconds":float(seconds)})
    ledgers=[line for line in log.splitlines() if "[RK4E]" in line]
    if not ledgers:
        raise ValueError("Missing native energy ledger")
    injection={}
    for key in ("floor_cum","de_cum","sfl_cum"):
        match=re.search(rf"\b{key}=([0-9.eE+-]+)",ledgers[-1])
        if match is None:
            raise ValueError(f"Missing {key} ledger")
        injection[key]=float(match.group(1))
    checks={
        "positive_state":bool(s0["positive"] and s1["positive"]),
        "fixed_mass_exact":bool(np.array_equal(p0["mass"],p1["mass"])),
        "initial_W0":bool(np.all(p0["w2"]==0)),
        "finite_nonnegative_new_weights":bool(np.all(weight>=0) and
                                               np.any(weight[active]>0)),
        "bounded_area_tiling":bool(abs(s1["volume_sum"]-AREA)<1e-5),
        "native_matches_independent_power_area":bool(np.max(area_error)<1e-4),
        "density_and_EOS":bool(max(s1["density_relative_max"],
                                    s1["eos_pressure_relative_max"])<1e-6),
        "short_total_energy":bool(energy_relative<1e-6),
        "short_incremental_orbit":bool(orbit_max<1e-4),
        "four_projected_RK4_stages":bool([x["stage"] for x in stages]==[0,1,2,3]
                                          and all(x["active"]==49038 for x in stages)),
        "zero_stage_weight_power":bool(len(stages)==4 and
                                        all(abs(x["projected_pressure_weight_power"])<1e-14
                                            for x in stages)),
        "weighted_variational_pressure_active":bool(any(marker in log for marker in (
            "[VARPRESS_W0] gated weighted fixed-mass variational pressure active",
            "[VARPRESS_W0] experimental smooth fixed-mass variational pressure active"))),
        "zero_injection":bool(all(x==0. for x in injection.values())),
        "no_geometry_abort":bool("[GPU_CLIP_CAPACITY]" not in log and
                                 "[GPU_TESS_INVALID]" not in log),
        "no_wall_crossing":bool("[FIXED_POLYGON_CROSSING]" not in log),
        "fixed_polygon_enabled":bool("[FIXED_POLYGON] enabled" in log),
    }
    result={
        "run":str(args.run.resolve()),"restart":str(restart.resolve()),
        "restart_sha256":digest(restart),"endpoint":str(endpoint.resolve()),
        "endpoint_sha256":digest(endpoint),"start_time":t0,"endpoint_time":args.time,
        "initial_state":s0,"final_state":s1,
        "weight_min":float(weight.min()),"weight_max":float(weight.max()),
        "active_weight_max":float(weight[active].max()),
        "independent_active_area_error_max":float(area_error.max()),
        "independent_active_area_error_p99":float(np.quantile(area_error,.99)),
        "energy_relative_change_from_restart":energy_relative,
        "disk_incremental_orbit_error_max":orbit_max,
        "rk4_injection":injection,"gpu_weight_stages":stages,
        "checks":checks,
        "scope":"One native fixed-mass weighted RK4 step from accepted t=6.5 W0 state. Global power and one-step geometry/energy audited; local material-face and long-time physics remain open.",
    }
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(result,indent=2,allow_nan=False))
    if not all(checks.values()):
        raise SystemExit("Native weighted one-step gate failed")


if __name__=="__main__":
    main()

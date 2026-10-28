"""Gate a native weighted bounded Kepler t=10 endpoint against its t=0 IC."""

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
DOMAIN_AREA = .5 * 1024 * 2.75**2 * np.sin(2*np.pi/1024)
APOTHEM = 2.75 * np.cos(np.pi/1024)
GAMMA_MINUS_ONE = float(np.float32(1.666667)) - 1.
STAGE_PATTERN = re.compile(r"\[GPU_VOLW2\] stage=(\d+) active=(\d+) faces=(\d+) "
                           r"alpha=([\deE+.-]+) work=([\deE+.-]+) elapsed=([\deE+.-]+) s")


def weighted_abs(value, mass):
    return float(np.sum(np.abs(value)*mass, dtype=np.longdouble)
                 / np.sum(mass, dtype=np.longdouble))


def energy(p):
    x = p["x"].astype(float)-3.
    y = p["y"].astype(float)-3.
    mass = p["mass"].astype(float)
    vx = p["vx"].astype(float)
    vy = p["vy"].astype(float)
    return float(np.sum(.5*mass*(vx*vx+vy*vy)+p["ie"].astype(float)
                        -mass/np.sqrt(x*x+y*y+.01175**2), dtype=np.longdouble))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("initial", type=Path)
    ap.add_argument("run", type=Path)
    ap.add_argument("--time", type=float, default=10.)
    ap.add_argument("--restart-time", type=float, default=6.501)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if abs(snapshot_time(args.initial,N))>1e-12:
        raise ValueError("Initial snapshot is not exact t=0")
    restart = args.run / "khout.000001.dat"
    if abs(snapshot_time(restart,N)-args.restart_time)>1e-10:
        raise ValueError("Restart is not at the requested accepted weighted time")
    endpoints = [path for path in args.run.glob("khout.*.dat")
                 if abs(snapshot_time(path,N)-args.time)<1e-10]
    if len(endpoints)!=1:
        raise ValueError("Expected unique complete exact endpoint")
    endpoint=endpoints[0]
    p0,_=ordered(args.initial,N)
    pr,_=ordered(restart,N)
    p,_=ordered(endpoint,N)
    s0,sr,s1=(state_metrics(q) for q in (p0,pr,p))
    x0=p0["x"].astype(float)-3.
    y0=p0["y"].astype(float)-3.
    r0=np.hypot(x0,y0)
    disk=(r0>=.5)&(r0<2.)
    ids=np.flatnonzero(disk)
    mass=p0["mass"][disk].astype(float)
    phase=args.time*(r0*r0+.01175**2)**(-.75)
    x_exact=x0*np.cos(phase)-y0*np.sin(phase)
    y_exact=x0*np.sin(phase)+y0*np.cos(phase)
    x=p["x"].astype(float)-3.
    y=p["y"].astype(float)-3.
    r=np.hypot(x,y)
    vx=p["vx"].astype(float)
    vy=p["vy"].astype(float)
    density_change=p["den"][disk].astype(float)/p0["den"][disk].astype(float)-1.
    entropy_log=(np.log(p["ie"][disk].astype(float)/p0["ie"][disk].astype(float))
                 +GAMMA_MINUS_ONE*np.log(p["volume"][disk].astype(float)
                                               /p0["volume"][disk].astype(float)))
    position_error=np.hypot(x[disk]-x_exact[disk],y[disk]-y_exact[disk])
    vr=(x*vx+y*vy)/np.maximum(r,1e-30)
    vc=r0*(r0*r0+.01175**2)**(-.75)
    radial_rms=float(np.sqrt(np.sum(mass*(vr[disk]/vc[disk])**2,dtype=np.longdouble)
                             /np.sum(mass,dtype=np.longdouble)))
    L0=float(np.sum(mass*(x0[disk]*p0["vy"][disk].astype(float)
                          -y0[disk]*p0["vx"][disk].astype(float)),dtype=np.longdouble))
    L1=float(np.sum(mass*(x[disk]*vy[disk]-y[disk]*vx[disk]),dtype=np.longdouble))
    weights=p["w2"].astype(float)
    lookup=np.full(N,-1,dtype=int)
    lookup[ids]=np.arange(len(ids))
    power_area,meta=power_areas(np.column_stack((x,y)),weights,ids,lookup)
    valid_area=np.isfinite(power_area)&(power_area>0)
    area_error=(np.abs(p["volume"][disk].astype(float)/power_area-1.)
                if np.all(valid_area) else None)
    area_error_max=float(np.max(area_error)) if area_error is not None else None
    e0,e1=energy(p0),energy(p)
    log=(args.run/"log").read_text(errors="replace")
    stages=[]
    for line in log.splitlines():
        match=STAGE_PATTERN.fullmatch(line)
        if match:
            st,nactive,faces,alpha,work,seconds=match.groups()
            stages.append({"stage":int(st),"active":int(nactive),"faces":int(faces),
                           "alpha":float(alpha),"work":float(work),"seconds":float(seconds)})
    ledgers=[line for line in log.splitlines() if "[RK4E]" in line]
    if not ledgers:
        raise ValueError("Missing native RK4 energy ledger")
    injection={}
    for key in ("floor_cum","de_cum","sfl_cum"):
        match=re.search(rf"\b{key}=([0-9.eE+-]+)",ledgers[-1])
        if not match:
            raise ValueError(f"Missing {key} energy ledger")
        injection[key]=float(match.group(1))
    density_abs=weighted_abs(density_change,mass)
    entropy_abs=weighted_abs(entropy_log,mass)
    orbit_abs=weighted_abs(position_error,mass)
    energy_rel=abs(e1-e0)/abs(e0)
    angular_drift=L1/L0-1.
    stage_sequence=bool(stages and len(stages)%4==0 and
                        all(row["stage"]==index%4 and row["active"]==49038
                            for index,row in enumerate(stages)))
    checks={
        "finite_positive":bool(s0["positive"] and sr["positive"] and s1["positive"]),
        "fixed_mass_exact":bool(np.array_equal(p0["mass"],pr["mass"]) and
                                np.array_equal(p0["mass"],p["mass"])),
        "initial_W0_and_weighted_restart":bool(np.all(p0["w2"]==0) and
                                                np.ptp(pr["w2"].astype(float))>1e-10),
        "finite_nonnegative_variable_weights":bool(np.all(np.isfinite(weights)) and
                                                     np.min(weights)>=0 and np.ptp(weights)>1e-10),
        "initial_and_final_area_tile_polygon":bool(abs(s0["volume_sum"]-DOMAIN_AREA)<1e-6 and
                                                    abs(s1["volume_sum"]-DOMAIN_AREA)<1e-5),
        "independent_power_area":bool(meta["hidden_active"]==0 and
                                      area_error_max is not None and area_error_max<1e-4),
        "density_and_EOS":bool(max(s0["density_relative_max"],s0["eos_pressure_relative_max"],
                                   sr["density_relative_max"],sr["eos_pressure_relative_max"],
                                   s1["density_relative_max"],s1["eos_pressure_relative_max"])<1e-6),
        "zero_injection":bool(all(value==0 for value in injection.values())),
        "total_energy":bool(energy_rel<1e-3),
        "disk_density":bool(density_abs<1e-3),
        "disk_entropy":bool(entropy_abs<1e-2),
        "disk_orbit":bool(orbit_abs<1e-2),
        "disk_radial_velocity":bool(radial_rms<.05),
        "disk_angular_momentum":bool(abs(angular_drift)<.01),
        "weight_stage_sequence":stage_sequence,
        "zero_stage_weight_power":bool(stage_sequence and
                                        all(abs(row["work"])<1e-14 for row in stages)),
        "inside_polygon_incircle":bool(np.all(r<APOTHEM)),
        "no_wall_crossing":bool("[FIXED_POLYGON_CROSSING]" not in log),
        "fixed_polygon_enabled":bool("[FIXED_POLYGON] enabled" in log),
    }
    result={
        "initial":str(args.initial.resolve()),"restart":str(restart.resolve()),
        "endpoint":str(endpoint.resolve()),"initial_sha256":digest(args.initial),
        "restart_sha256":digest(restart),"endpoint_sha256":digest(endpoint),
        "time":args.time,"restart_time":args.restart_time,
        "n":N,"disk_sites":int(disk.sum()),
        "initial_state":s0,"restart_state":sr,"endpoint_state":s1,
        "total_energy_relative":energy_rel,"rk4_injection":injection,
        "disk_density_abs_mass_weighted":density_abs,
        "disk_density_p99":float(np.quantile(np.abs(density_change),.99)),
        "disk_entropy_log_abs_mass_weighted":entropy_abs,
        "disk_position_error_abs_mass_weighted":orbit_abs,
        "disk_radial_speed_rms_over_circular":radial_rms,
        "disk_angular_momentum_relative":angular_drift,
        "weight_min":float(np.min(weights)),"weight_max":float(np.max(weights)),
        "weight_range":float(np.ptp(weights)),
        "independent_power_area_error_max":area_error_max,
        "independent_power_area_error_p99":(float(np.quantile(area_error,.99))
                                             if area_error is not None else None),
        "hidden_active_power_cells":meta["hidden_active"],
        "weight_stages":len(stages),
        "max_abs_stage_weight_power":max((abs(row["work"]) for row in stages),default=None),
        "mean_stage_rate_seconds":float(np.mean([row["seconds"] for row in stages])) if stages else None,
        "max_current_site_radius":float(np.max(r)),
        "checks":checks,
        "scope":"Experimental bounded one-rank fixed-mass weighted Kepler endpoint. The exact-orbit density area reference for W0 is intentionally not used for weighted cells; GIZMO t=600 remains open.",
    }
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(result,indent=2,allow_nan=False))
    if not all(checks.values()):
        raise SystemExit("Weighted native Kepler t=10 physics gate failed")


if __name__=="__main__":
    main()

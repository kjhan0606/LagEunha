"""Gate one native weighted RK4 step started at the accepted jitter t=0 IC."""

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
GAMMA_MINUS_ONE = float(np.float32(1.666667))-1.
DOMAIN_AREA = .5*1024*2.75**2*np.sin(2*np.pi/1024)
APOTHEM = 2.75*np.cos(np.pi/1024)
FIELDS = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2")


def energy(p):
    x=p["x"].astype(float)-3.
    y=p["y"].astype(float)-3.
    mass=p["mass"].astype(float)
    vx=p["vx"].astype(float)
    vy=p["vy"].astype(float)
    return float(np.sum(.5*mass*(vx*vx+vy*vy)+p["ie"].astype(float)
                        -mass/np.sqrt(x*x+y*y+.01175**2), dtype=np.longdouble))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("run", type=Path)
    ap.add_argument("--expected-initial-sha256", required=True)
    ap.add_argument("--expected-binary-sha256", required=True)
    ap.add_argument("--output", type=Path, required=True)
    args=ap.parse_args()
    reference=args.run/"accepted_initial.dat"
    start=args.run/"khout.000000.dat"
    endpoint=args.run/"khout.000001.dat"
    if digest(reference)!=args.expected_initial_sha256:
        raise ValueError("Accepted initial reference hash changed")
    binary_sha=digest(args.run/"eunha2.exe")
    if binary_sha!=args.expected_binary_sha256:
        raise ValueError("Weighted executable hash changed")
    if abs(snapshot_time(start,N))>1e-12 or abs(snapshot_time(endpoint,N)-.001)>1e-10:
        raise ValueError("Expected complete native t=0 and t=.001 snapshots")
    ref,_=ordered(reference,N)
    p0,_=ordered(start,N)
    p1,_=ordered(endpoint,N)
    s0,s1=state_metrics(p0),state_metrics(p1)
    same_fields={field:bool(np.array_equal(ref[field],p0[field])) for field in FIELDS}
    xy0=np.column_stack((p0["x"].astype(float)-3.,p0["y"].astype(float)-3.))
    r0=np.linalg.norm(xy0,axis=1)
    disk=(r0>=.5)&(r0<2.)
    ids=np.flatnonzero(disk)
    mass=p0["mass"][disk].astype(float)
    xy1=np.column_stack((p1["x"].astype(float)-3.,p1["y"].astype(float)-3.))
    weight=p1["w2"].astype(float)
    lookup=np.full(N,-1,dtype=int)
    lookup[ids]=np.arange(len(ids))
    area,meta=power_areas(xy1,weight,ids,lookup)
    area_error=(float(np.max(np.abs(p1["volume"][ids].astype(float)/area-1.)))
                if meta["hidden_active"]==0 and np.all(area>0) else float("inf"))
    drift=p1["den"][disk].astype(float)/p0["den"][disk].astype(float)-1.
    density_abs=float(np.average(np.abs(drift),weights=mass))
    entropy=np.log(p1["ie"][disk].astype(float)/p0["ie"][disk].astype(float))
    entropy+=GAMMA_MINUS_ONE*np.log(p1["volume"][disk].astype(float)
                                      /p0["volume"][disk].astype(float))
    entropy_abs=float(np.average(np.abs(entropy),weights=mass))
    e0,e1=energy(p0),energy(p1)
    log=(args.run/"log").read_text(errors="replace")
    stages=[int(m.group(1)) for m in re.finditer(r"^\[GPU_VOLW2\] stage=(\d+)",log,re.M)]
    ledgers=[line for line in log.splitlines() if line.startswith("[RK4E]")]
    injection={}
    if ledgers:
        for name in ("floor_cum","de_cum","sfl_cum"):
            found=re.search(rf"\b{name}=([0-9.eE+-]+)",ledgers[-1])
            injection[name]=float(found.group(1)) if found else float("nan")
    checks={
        "accepted_t0_physical_fields_exact":all(same_fields.values()),
        "finite_positive":bool(s0["positive"] and s1["positive"]),
        "fixed_mass_exact":bool(np.array_equal(p0["mass"],p1["mass"])),
        "initial_W0":bool(np.all(p0["w2"]==0)),
        "nonnegative_nontrivial_weights":bool(np.all(np.isfinite(weight)) and
                                              np.min(weight)>=0 and np.ptp(weight)>1e-12),
        "area_tiles_polygon":bool(abs(s0["volume_sum"]-DOMAIN_AREA)<1e-6 and
                                   abs(s1["volume_sum"]-DOMAIN_AREA)<1e-5),
        "independent_power_area":bool(area_error<1e-4),
        "density_and_EOS":bool(max(s0["density_relative_max"],
                                    s0["eos_pressure_relative_max"],
                                    s1["density_relative_max"],
                                    s1["eos_pressure_relative_max"])<1e-6),
        "one_step_density_safety":bool(density_abs<.01),
        "one_step_entropy_safety":bool(entropy_abs<.01),
        "total_energy":bool(abs(e1-e0)/abs(e0)<1e-6),
        "zero_injection":bool(len(injection)==3 and all(x==0 for x in injection.values())),
        "four_weight_stages":stages==[0,1,2,3],
        "pressure_force_active":bool("[VARPRESS_W0]" in log),
        "inside_polygon_incircle":bool(np.all(np.linalg.norm(xy1,axis=1)<APOTHEM)),
        "no_geometry_abort":bool("[GPU_TESS_INVALID]" not in log and
                                  "[GPU_CLIP_CAPACITY]" not in log),
        "no_wall_crossing":bool("[FIXED_POLYGON_CROSSING]" not in log),
        "fixed_polygon_enabled":bool("[FIXED_POLYGON] enabled" in log),
    }
    result={
        "reference":str(reference.resolve()),"initial":str(start.resolve()),
        "endpoint":str(endpoint.resolve()),"reference_sha256":digest(reference),
        "initial_sha256":digest(start),"endpoint_sha256":digest(endpoint),
        "binary_sha256":binary_sha,"time":.001,"n":N,"disk_sites":int(disk.sum()),
        "initial_field_matches":same_fields,"initial_state":s0,"endpoint_state":s1,
        "disk_density_abs_mass_weighted":density_abs,
        "disk_entropy_log_abs_mass_weighted":entropy_abs,
        "total_energy_relative":abs(e1-e0)/abs(e0),
        "independent_power_area_error_max":area_error,
        "hidden_active_power_cells":meta["hidden_active"],
        "weight_range":float(np.ptp(weight)),"weight_stages":stages,
        "rk4_injection":injection,"checks":checks,
        "scope":"Native one-step fresh weighted jitter pilot. Passing is not a "
                "long-time or local moving-face conservation result.",
    }
    with args.output.open("x") as output:
        json.dump(result,output,indent=2,allow_nan=False)
        output.write("\n")
    print(json.dumps({"failed":[k for k,v in checks.items() if not v],
                      "density_abs":density_abs,"area_error_max":area_error},indent=2))
    if not all(checks.values()):
        raise SystemExit("Fresh weighted t=.001 native pilot failed")


if __name__=="__main__":
    main()

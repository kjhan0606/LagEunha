"""One-step causal-control comparison; explicitly not a physical hydro gate."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


N=55351
BASE=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/force_inputs_t65_a10_v1/jitter")
DT=.0005


def rms(vector,mass):
    return float(np.sqrt(np.sum(mass*vector*vector,dtype=np.longdouble)
                         /np.sum(mass,dtype=np.longdouble)))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("run",type=Path)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    run=args.run.resolve()
    staged=json.loads((run/"staging.json").read_text())
    original=json.loads((BASE/"audit_short.json").read_text())
    diagnostic=json.loads((run/"audit_short.json").read_text())
    original_end=Path(original["endpoint"])
    diagnostic_end=Path(diagnostic["endpoint"])
    if (not all(original["checks"].values()) or not all(diagnostic["checks"].values())
            or original["restart_sha256"]!=diagnostic["restart_sha256"]
            or original["restart_sha256"]!=staged["accepted_parent_endpoint_sha256"]
            or digest(original_end)!=original["endpoint_sha256"]
            or digest(diagnostic_end)!=diagnostic["endpoint_sha256"]
            or digest(run/"eunha2.exe")!=staged["diagnostic_binary_sha256"]):
        raise ValueError("Changed or failed one-step native inputs/outputs")
    for name in ("params.dat","glass.bin","areas.bin"):
        if digest(BASE/name)!=digest(run/name):
            raise ValueError(f"Physical input differs from original: {name}")
    log=(run/"log").read_text(errors="replace")
    if "[DIAG_PRESSURE_FORCE_OFF] pressure force suppressed" not in log:
        raise ValueError("Diagnostic pressure-force-off marker missing")
    p0,_=ordered(run/"khout.000001.dat",N)
    po,mo=ordered(original_end,N)
    pd,md=ordered(diagnostic_end,N)
    if abs(float(mo["t"])-6.5005)>1e-10 or abs(float(md["t"])-6.5005)>1e-10:
        raise ValueError("One-step controls have different endpoint time")
    r0=np.hypot(p0["x"].astype(float)-3.,p0["y"].astype(float)-3.)
    disk=(r0>=.5)&(r0<2.)
    mass=p0["mass"][disk].astype(float)
    velocity_difference=np.hypot(po["vx"][disk]-pd["vx"][disk],
                                 po["vy"][disk]-pd["vy"][disk])
    position_difference=np.hypot(po["x"][disk]-pd["x"][disk],
                                 po["y"][disk]-pd["y"][disk])
    density_difference=np.abs(po["den"][disk].astype(float)-pd["den"][disk].astype(float))
    force_reference=json.loads((run/"pressure_force_t65.json").read_text())
    if force_reference["endpoint_sha256"] != staged["accepted_parent_endpoint_sha256"]:
        raise ValueError("Frozen force reference uses a different t6.5 checkpoint")
    force_rms=force_reference["whole_disk"]["acceleration_mass_weighted_rms"]
    velocity_rms=rms(velocity_difference,mass)
    result={"original_one_step_audit":str(BASE/"audit_short.json"),
            "diagnostic_one_step_audit":str(run/"audit_short.json"),
            "accepted_parent_sha256":staged["accepted_parent_endpoint_sha256"],
            "original_endpoint_sha256":original["endpoint_sha256"],
            "diagnostic_endpoint_sha256":diagnostic["endpoint_sha256"],
            "diagnostic_binary_sha256":staged["diagnostic_binary_sha256"],
            "time":6.5005,"dt":DT,
            "disk_velocity_difference_mass_weighted_rms":velocity_rms,
            "disk_velocity_difference_over_dt_rms":velocity_rms/DT,
            "frozen_t65_pressure_force_rms":force_rms,
            "difference_acceleration_over_frozen_pressure_force_rms":velocity_rms/DT/force_rms,
            "disk_position_difference_mass_weighted_rms":rms(position_difference,mass),
            "disk_density_difference_mass_weighted_mean":float(
                np.sum(mass*density_difference,dtype=np.longdouble)/np.sum(mass,dtype=np.longdouble)),
            "original_short_energy_relative":original["energy_relative_change_from_restart"],
            "diagnostic_short_energy_relative":diagnostic["energy_relative_change_from_restart"],
            "original_short_orbit_error_max":original["disk_incremental_orbit_error_max"],
            "diagnostic_short_orbit_error_max":diagnostic["disk_incremental_orbit_error_max"],
            "checks":{"same_restart_sha256":True,"both_short_safety_gates":True,
                      "changed_disk_velocity":bool(velocity_rms>1e-12),
                      "bounded_intervention":bool(velocity_rms/DT<1e-4)},
            "scope":"Pressure-force-off causal control at one saved RK4 step. Geometric PdV/EOS retained; energy conservation is not the physical variational operator. No long-time repair claim."}
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(result,indent=2))
    if not all(result["checks"].values()):
        raise SystemExit("One-step pressure-force-off intervention gate failed")


if __name__=="__main__":
    main()

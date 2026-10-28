"""Compare native raw/projected moving weights from the same accepted t=0 IC."""

import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_q2_weighted_fresh_pilot_20261027 import energy
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DOCS = Path(__file__).resolve().parent
N = 55351
MODES = ("projected", "raw")
OUTPUT = DOCS / "KEPLER_Q2_FRESH_PROJECTION_PAIR_T001_2026-10-27.json"


def ring_metrics(initial, endpoint, ring, order, disk):
    mass = initial["mass"].astype(float)
    x0, y0 = initial["x"].astype(float)-3, initial["y"].astype(float)-3
    x1, y1 = endpoint["x"].astype(float)-3, endpoint["y"].astype(float)-3
    l0 = x0*initial["vy"].astype(float)-y0*initial["vx"].astype(float)
    l1 = x1*endpoint["vy"].astype(float)-y1*endpoint["vx"].astype(float)
    pressure = endpoint["pressure"].astype(float)
    p0 = initial["pressure"].astype(float)
    q = pressure[order]/p0[order]
    angular_jump = q-np.roll(q,-1)
    volume = endpoint["volume"].astype(float)/initial["volume"].astype(float)-1
    density = endpoint["den"].astype(float)/initial["den"].astype(float)-1
    return {
        "first_ring_count": int(ring.sum()),
        "first_ring_volume_change_rms": float(np.sqrt(np.mean(volume[ring]**2))),
        "first_ring_density_change_rms": float(np.sqrt(np.mean(density[ring]**2))),
        "first_ring_pressure_change_mean_abs": float(np.mean(np.abs(pressure[ring]/p0[ring]-1))),
        "first_ring_angular_pressure_jump_rms": float(np.sqrt(np.mean(angular_jump**2))),
        "first_ring_orbital_ell_change_mean_abs": float(np.average(
            np.abs(l1[ring]/l0[ring]-1), weights=mass[ring])),
        "disk_density_change_mean_abs": float(np.average(np.abs(density[disk]),weights=mass[disk])),
        "total_energy_relative_change": float((energy(endpoint)-energy(initial))/abs(energy(initial))),
    }


def main():
    result = {"time": .001, "modes": {}, "scope":
              "Same-binary, same-fresh-t0 native causal one-step comparison; "
              "not a long-time Kepler stability claim."}
    initial_ref = None
    binary_ref = None
    for mode in MODES:
        run = ROOT / f"gpu_energy_neutral_v11b_fresh_{mode}_t001_a10_v1"
        manifest = json.loads((run/"stage_manifest.json").read_text())
        audit = json.loads((run/"audit_t001.json").read_text())
        if (manifest["mode"] != mode or manifest["target_time"] != .001
                or len(audit["checks"]) != 18 or not all(audit["checks"].values())
                or audit["binary_sha256"] != manifest["binary_sha256"]):
            raise ValueError(f"Unaccepted native {mode} endpoint")
        start, end = run/"khout.000000.dat", run/"khout.000001.dat"
        if (digest(start) != audit["initial_sha256"]
                or digest(end) != audit["endpoint_sha256"]
                or abs(snapshot_time(start,N)) > 1e-12
                or abs(snapshot_time(end,N)-.001) > 1e-10):
            raise ValueError(f"Changed or incomplete native {mode} snapshot")
        initial,_ = ordered(start,N)
        endpoint,_ = ordered(end,N)
        if initial_ref is None:
            initial_ref = initial
            binary_ref = audit["binary_sha256"]
        elif (audit["binary_sha256"] != binary_ref or
              any(not np.array_equal(initial[name],initial_ref[name])
                  for name in initial.dtype.names)):
            raise ValueError("The pair does not share the same initial state and binary")
        x0,y0=initial["x"].astype(float)-3,initial["y"].astype(float)-3
        r0=np.hypot(x0,y0)
        disk=(r0>=.5)&(r0<2)
        ring=(r0>=.5)&(r0<.52)
        ring_ids=np.flatnonzero(ring)
        order=ring_ids[np.argsort(np.arctan2(y0[ring_ids],x0[ring_ids]))]
        log=(run/"log").read_text(errors="replace")
        works=[float(x) for x in re.findall(r"^\[GPU_VOLW2\] stage=\d+.*?\bwork=([\deE+.-]+)",log,re.M)]
        if len(works) != 4:
            raise ValueError(f"Expected four RK4 weight works for {mode}, got {len(works)}")
        result["modes"][mode]={
            "run":str(run),"initial_sha256":audit["initial_sha256"],
            "endpoint_sha256":audit["endpoint_sha256"],
            "binary_sha256":audit["binary_sha256"],
            "audit_passed":len(audit["checks"]),
            "gpu_stage_pressure_weight_works":works,
            **ring_metrics(initial,endpoint,ring,order,disk),
        }
    p=result["modes"]["projected"]
    r=result["modes"]["raw"]
    result["initial_state_fields_exact"] = list(initial_ref.dtype.names)
    result["initial_file_hashes_differ"] = p["initial_sha256"] != r["initial_sha256"]
    result["ratios_projected_over_raw"]={
        key:p[key]/r[key] for key in (
            "first_ring_volume_change_rms", "first_ring_density_change_rms",
            "first_ring_pressure_change_mean_abs", "first_ring_angular_pressure_jump_rms",
            "first_ring_orbital_ell_change_mean_abs") if r[key] != 0
    }
    with OUTPUT.open("x") as output:
        json.dump(result,output,indent=2,allow_nan=False)
        output.write("\n")
    print(json.dumps({"output":str(OUTPUT),"ratios":result["ratios_projected_over_raw"]},indent=2))


if __name__ == "__main__":
    main()

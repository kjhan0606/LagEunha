"""Estimate local pressure-weight work in the accepted native weighted step.

The estimate uses midpoint snapshot states, not every native RK4 stage; its
purpose is to reveal local work hidden by a zero global projection.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_weight_work_t65_20261027 import faces_touching_active
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
DT = .001


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=("aligned", "jitter"))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    base = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
    run = base / "gpu_energy_neutral_weighted_t6501_guard_a10_v2" / args.case
    receipt = json.loads((run / "audit_t6501.json").read_text())
    if not all(receipt["checks"].values()):
        raise ValueError("Weighted step did not pass its native gate")
    start_path = Path(receipt["restart"])
    end_path = Path(receipt["endpoint"])
    if (digest(start_path) != receipt["restart_sha256"]
            or digest(end_path) != receipt["endpoint_sha256"]):
        raise ValueError("Changed native snapshots")
    start, _ = ordered(start_path, N)
    end, _ = ordered(end_path, N)
    xy0 = np.column_stack((start["x"].astype(float)-3., start["y"].astype(float)-3.))
    xy1 = np.column_stack((end["x"].astype(float)-3., end["y"].astype(float)-3.))
    xy = .5*(xy0+xy1)
    weight0 = start["w2"].astype(float)
    weight1 = end["w2"].astype(float)
    weight = .5*(weight0+weight1)
    rate = (weight1-weight0)/DT
    pressure = .5*(start["pressure"].astype(float)+end["pressure"].astype(float))
    internal = .5*(start["ie"].astype(float)+end["ie"].astype(float))
    r0 = np.linalg.norm(xy0, axis=1)
    active = np.flatnonzero((r0>=.5)&(r0<2.))
    active_flag = np.zeros(N, dtype=bool)
    active_flag[active] = True
    lookup = np.full(N, -1, dtype=int)
    lookup[active] = np.arange(len(active))
    area, meta, triangles, centers = power_areas(
        xy, weight, active, lookup, return_geometry=True)
    if meta["hidden_active"] or not np.all(area>0):
        raise ValueError("Hidden midpoint power cell")
    i,j,length = faces_touching_active(triangles, centers, active_flag)
    distance = np.linalg.norm(xy[j]-xy[i], axis=1)
    if not np.all(distance>0):
        raise ValueError("Coincident midpoint sites")
    g = length/(2*distance)
    dv = g*(rate[i]-rate[j])
    cell_dv = np.zeros(N)
    np.add.at(cell_dv, i, dv)
    np.add.at(cell_dv, j, -dv)
    local_power = pressure*cell_dv
    face_power = (pressure[i]-pressure[j])*dv
    cell_sum = float(np.sum(local_power, dtype=np.longdouble))
    face_sum = float(np.sum(face_power, dtype=np.longdouble))
    face_abs = float(np.sum(np.abs(face_power), dtype=np.longdouble))
    if abs(cell_sum-face_sum)>1e-12*max(face_abs,1e-30):
        raise ValueError("Face/cell pressure-work identity failed")
    frac = DT*np.abs(local_power[active])/internal[active]
    result = {
        "case": args.case, "start_sha256": receipt["restart_sha256"],
        "end_sha256": receipt["endpoint_sha256"], "dt": DT,
        "midpoint_disk_area_error_max": float(np.max(np.abs(
            .5*(start["volume"][active].astype(float)+end["volume"][active].astype(float))/area-1))),
        "disk_touching_faces": len(i),
        "signed_pressure_weight_power_rate": face_sum,
        "absolute_face_pressure_weight_power_rate": face_abs,
        "face_cell_identity_error": abs(cell_sum-face_sum),
        "one_step_signed_work_over_domain_internal_energy":
            DT*face_sum/float(np.sum(internal, dtype=np.longdouble)),
        "one_step_absolute_face_work_over_domain_internal_energy":
            DT*face_abs/float(np.sum(internal, dtype=np.longdouble)),
        "active_local_weight_work_over_internal_energy_p50": float(np.quantile(frac,.5)),
        "active_local_weight_work_over_internal_energy_p99": float(np.quantile(frac,.99)),
        "active_local_weight_work_over_internal_energy_max": float(np.max(frac)),
        "scope": "Midpoint snapshot estimate of local pressure-weight work. "
                 "Global RK4 projection is separately recorded in the native stage log; "
                 "this estimate is not an exact stage-by-stage local first-law audit.",
    }
    with args.output.open("x") as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

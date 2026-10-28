"""Estimate material-face normal slip of the actual native weighted step."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_weighted_face_slip_t65_20261027 import diagram
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
DT = .001
INITIAL_RECEIPTS = {
    "aligned": "KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json",
    "jitter": "KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=("aligned", "jitter"))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    docs = Path(__file__).resolve().parent
    base = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
    run = base / "gpu_energy_neutral_weighted_t6501_guard_a10_v2" / args.case
    receipt = json.loads((run / "audit_t6501.json").read_text())
    initial_receipt = json.loads((docs / INITIAL_RECEIPTS[args.case]).read_text())
    if (not all(receipt["checks"].values())
            or receipt["restart_sha256"] != initial_receipt["endpoint_sha256"]):
        raise ValueError("Unaccepted native lineage")
    for path, expected in ((Path(receipt["restart"]), receipt["restart_sha256"]),
                           (Path(receipt["endpoint"]), receipt["endpoint_sha256"]),
                           (Path(initial_receipt["initial"]), initial_receipt["initial_sha256"])):
        if digest(path) != expected:
            raise ValueError(f"Changed snapshot: {path}")
    initial, _ = ordered(Path(initial_receipt["initial"]), N)
    start, _ = ordered(Path(receipt["restart"]), N)
    end, _ = ordered(Path(receipt["endpoint"]), N)
    xy0 = np.column_stack((initial["x"].astype(float)-3., initial["y"].astype(float)-3.))
    xy_start = np.column_stack((start["x"].astype(float)-3., start["y"].astype(float)-3.))
    xy_end = np.column_stack((end["x"].astype(float)-3., end["y"].astype(float)-3.))
    xy = .5*(xy_start+xy_end)
    velocity = .5*(
        np.column_stack((start["vx"].astype(float), start["vy"].astype(float)))+
        np.column_stack((end["vx"].astype(float), end["vy"].astype(float))))
    r0 = np.linalg.norm(xy0, axis=1)
    active = np.flatnonzero((r0>=.5)&(r0<2.))
    lookup = np.full(N,-1,dtype=int)
    lookup[active] = np.arange(len(active))
    weight0 = start["w2"].astype(float)
    weight1 = end["w2"].astype(float)
    weight = .5*(weight0+weight1)
    rate = (weight1-weight0)/DT
    zero = np.zeros(N)
    baseline = diagram(xy, velocity, r0, active, lookup, zero, zero)
    weighted_fixed = diagram(xy, velocity, r0, active, lookup, weight, zero)
    weighted_moving = diagram(xy, velocity, r0, active, lookup, weight, rate)
    reference = baseline["all_disk_disk"]["fixed"]["length_weighted_rms_slip"]
    actual = weighted_moving["all_disk_disk"]["moving"]["length_weighted_rms_slip"]
    result = {
        "case": args.case,
        "initial_sha256": initial_receipt["initial_sha256"],
        "start_sha256": receipt["restart_sha256"],
        "end_sha256": receipt["endpoint_sha256"],
        "midpoint_native_step_estimate": True,
        "W0_on_same_midpoint_sites": baseline,
        "weighted_geometry_frozen_W": weighted_fixed,
        "weighted_geometry_applied_rate": weighted_moving,
        "all_disk_disk_rms_slip_ratio_weighted_over_W0": actual/reference,
        "scope": "Analytic Kepler fluid normal speed versus midpoint geometric face speed. "
                 "The rate is an endpoint difference, not each exact native RK4 stage rate; "
                 "this does not prove material faces or native mass-flux conservation.",
    }
    with args.output.open("x") as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write("\n")
    print(json.dumps({"case": args.case,
                      "W0_rms": reference,
                      "weighted_rms": actual,
                      "ratio": actual/reference}, indent=2))


if __name__ == "__main__":
    main()

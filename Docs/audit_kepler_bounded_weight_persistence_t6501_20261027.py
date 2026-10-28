"""Frozen t=6.5 static Laguerre weights on accepted t=6.501 positions."""

import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_weight_feasibility_t65_20261027 import guard, metrics
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
STATIC = ROOT/"weight_feasibility_t65_cpu_v1"
PROBE = ROOT/"faceforce_pair_t6501_a10_v1"
CASES = {
    "aligned": DOCS/"KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json",
    "jitter": DOCS/"KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json",
}
GAMMA = float(np.float32(1.666667))


def mean_abs(value, mass):
    return float(np.sum(np.abs(value)*mass, dtype=np.longdouble)
                 / np.sum(mass, dtype=np.longdouble))


def one_case(label, parent_path):
    parent = json.loads(parent_path.read_text())
    short = json.loads((PROBE/label/"audit_t6501.json").read_text())
    static_path = STATIC/f"{label}_weight.json"
    static = json.loads(static_path.read_text())
    weights_path = STATIC/f"{label}_weight.npz"
    if (not all(parent["checks"].values()) or not all(short["checks"].values())
            or not static["converged_to_1e_minus_7"]
            or abs(parent["time"]-6.5)>1e-10
            or abs(short["endpoint_time"]-6.501)>1e-10
            or static["endpoint_sha256"] != parent["endpoint_sha256"]
            or short["restart_sha256"] != parent["endpoint_sha256"]):
        raise ValueError(f"{label}: accepted parent/static/short receipts disagree")
    for path, expected in ((Path(parent["initial"]),parent["initial_sha256"]),
                           (Path(parent["endpoint"]),parent["endpoint_sha256"]),
                           (Path(short["endpoint"]),short["endpoint_sha256"])):
        if digest(path) != expected:
            raise ValueError(f"{label}: snapshot differs from accepted SHA")
    p0,_ = ordered(Path(parent["initial"]),55351)
    p65,_ = ordered(Path(parent["endpoint"]),55351)
    p6501,_ = ordered(Path(short["endpoint"]),55351)
    with np.load(weights_path) as archive:
        active = archive["active"].astype(int)
        active_weight = archive["weight"].astype(float)
    radius0=np.hypot(p0["x"].astype(float)-3.,p0["y"].astype(float)-3.)
    expected_active=np.flatnonzero((radius0>=.5)&(radius0<2.))
    if not np.array_equal(active,expected_active) or len(active)!=49038:
        raise ValueError(f"{label}: changed active weight IDs")
    weight=np.zeros(55351)
    weight[active]=active_weight
    lookup=np.full(55351,-1,dtype=int)
    lookup[active]=np.arange(len(active))
    target=p0["volume"][active].astype(float)
    mass=p0["mass"][active].astype(float)
    result={
        "accepted_parent_sha256":parent["endpoint_sha256"],
        "accepted_one_step_sha256":short["endpoint_sha256"],
        "static_weight_json_sha256":digest(static_path),
        "static_weight_npz_sha256":digest(weights_path),
        "static_guard_max_ratio":static["native_face_guard"]["max_absolute_weight_difference_over_d2"],
    }
    for time,p in (("t65",p65),("t6501",p6501)):
        xy=np.column_stack((p["x"].astype(float)-3.,p["y"].astype(float)-3.))
        area,geom,triangles,_=power_areas(xy,weight,active,lookup,return_geometry=True)
        if geom["hidden_active"] or np.any(~np.isfinite(area)) or np.any(area<=0):
            raise ValueError(f"{label}: hidden or invalid frozen-weight cell at {time}")
        weighted_error=metrics(area,target,radius0[active])
        unweighted_area=p["volume"][active].astype(float)
        unweighted_error=metrics(unweighted_area,target,radius0[active])
        guard_state=guard(xy,weight,triangles,lookup)
        result[time]={
            "weighted_target_error":weighted_error,
            "unweighted_target_error":unweighted_error,
            "weighted_pressure_abs_relative_mass_weighted":mean_abs((target/area)**GAMMA-1.,mass),
            "unweighted_pressure_abs_relative_mass_weighted":mean_abs((target/unweighted_area)**GAMMA-1.,mass),
            "face_guard_max_ratio":guard_state["max_absolute_weight_difference_over_d2"],
            "faces_guard_over_0p8":guard_state["faces_over_0p8"],
            "faces_guard_over_1":guard_state["faces_over_1"],
            "hidden_active":geom["hidden_active"],
        }
    if abs(result["t65"]["weighted_target_error"]["max_relative"]-
           static["final_target_error"]["max_relative"]) > 1e-8:
        raise ValueError(f"{label}: static solution not reproduced")
    return result


def main():
    cases={label:one_case(label,path) for label,path in CASES.items()}
    result={
        "time_start":6.5,"time_next":6.501,
        "cases":cases,
        "scope":("Geometry-only reuse of t=6.5 static exact-volume weights on native W0 "
                 "one-step positions. The gas did not evolve with these weights; no "
                 "weighted face speed, force, work, mass or GPU runtime is validated."),
    }
    out=DOCS/"KEPLER_Q2_BOUNDED_WEIGHT_PERSIST_T6501_2026-10-27.json"
    with out.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    for label,q in cases.items():
        print(label,"t6501_weighted_p99",q["t6501"]["weighted_target_error"]["p99_relative"],
              "t6501_unweighted_p99",q["t6501"]["unweighted_target_error"]["p99_relative"],
              "weighted_pressure_abs",q["t6501"]["weighted_pressure_abs_relative_mass_weighted"],
              "unweighted_pressure_abs",q["t6501"]["unweighted_pressure_abs_relative_mass_weighted"],
              "guard_max",q["t6501"]["face_guard_max_ratio"])


if __name__=="__main__":
    main()

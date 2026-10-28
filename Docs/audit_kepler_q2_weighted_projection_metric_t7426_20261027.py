"""Compare energy-neutral weight-rate projections on one native power diagram."""

import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator, cg

from audit_kepler_sparse_area_t0825_20261027 import ordered
from export_kepler_q2_weight_rate_csr_t65_20261027 import finite_power_edges
from laguerre_newton_volume_t6_20261026 import newton_jacobian
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from laguerre_weight_rate_actual_t625_20261026 import face_volume_rate
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_weighted_t10_a10_v1")
OUTPUT = DOCS / "KEPLER_Q2_WEIGHTED_PROJECTION_METRIC_T7426_2026-10-27.json"
N = 55351
DT = .00193983


def jacobi(matrix, diag, rhs, count):
    rate = np.zeros(len(rhs))
    for _ in range(count):
        rate += .8*(rhs-matrix@rate)/diag
    return rate


def score(matrix, fixed, h, area, inner, initial, direction):
    denom = float(h@direction)
    if abs(denom) < 1e-35:
        raise ValueError("Nearly singular projection direction")
    alpha = float(h@initial)/denom
    projected = initial-alpha*direction
    residual = fixed+matrix@projected
    fractional = DT*np.abs(residual)/area
    return {
        "alpha": alpha,
        "energy_neutral_work": float(h@projected),
        "relative_rate_change_l2": float(np.linalg.norm(projected-initial)/np.linalg.norm(initial)),
        "volume_rate_relative_residual_l2": float(np.linalg.norm(residual)/np.linalg.norm(fixed)),
        "volume_change_disk_p99": float(np.quantile(fractional,.99)),
        "volume_change_inner_p99": float(np.quantile(fractional[inner],.99)),
        "volume_change_disk_mean_abs": float(np.mean(fractional)),
        "rate_abs_max": float(np.max(np.abs(projected))),
    }


def case(label):
    started = time.monotonic()
    run = ROOT / label
    initial_file = run / "initial.dat"
    checkpoint = run / "khout.000478.dat"
    t = snapshot_time(checkpoint,N)
    if not 7.42 < t < 7.43:
        raise ValueError("Wrong native checkpoint")
    p0,_ = ordered(initial_file,N)
    p,_ = ordered(checkpoint,N)
    xy0 = np.column_stack((p0["x"].astype(float)-3.,p0["y"].astype(float)-3.))
    xy = np.column_stack((p["x"].astype(float)-3.,p["y"].astype(float)-3.))
    r0 = np.linalg.norm(xy0,axis=1)
    active = np.flatnonzero((r0>=.5)&(r0<2.))
    inner = r0[active]<.75
    lookup = np.full(N,-1,dtype=int)
    lookup[active] = np.arange(len(active))
    W = p["w2"].astype(float)
    area,meta,triangles,centers = power_areas(xy,W,active,lookup,return_geometry=True)
    if meta["hidden_active"] or np.any(area<=0):
        raise ValueError("Invalid native power geometry")
    velocity = np.column_stack((p["vx"].astype(float),p["vy"].astype(float)))
    fixed,_,nfaces = face_volume_rate(xy,velocity,W,np.zeros(N),active,lookup,
                                       triangles,centers)
    matrix,diag = newton_jacobian(xy,lookup,len(active),triangles,centers)
    i,j,c1,c2 = finite_power_edges(triangles,centers,lookup>=0)
    if len(i)!=nfaces:
        raise ValueError("Power face graph mismatch")
    g = np.linalg.norm(c2-c1,axis=1)/(2*np.linalg.norm(xy[i]-xy[j],axis=1))
    pressure = p["pressure"].astype(float)
    full_h = np.zeros(N)
    pair = g*(pressure[i]-pressure[j])
    np.add.at(full_h,i,pair)
    np.add.at(full_h,j,-pair)
    h = full_h[active]
    if not np.linalg.norm(h)>0:
        raise ValueError("No pressure weight gradient")
    rhs=-fixed
    d16 = jacobi(matrix,diag,h,16)
    d64 = jacobi(matrix,diag,h,64)
    preconditioner=LinearOperator(matrix.shape,matvec=lambda z:z/diag)
    iterations=[0]
    d_exact,info=cg(matrix,h,M=preconditioner,rtol=1e-9,atol=0.,maxiter=2000,
                    callback=lambda _:iterations.__setitem__(0,iterations[0]+1))
    if info:
        raise RuntimeError(f"Pressure projection solve failed: {info}")
    directions={"euclidean_h":h,"jacobi16_J_inverse_h":d16,
                "jacobi64_J_inverse_h":d64,"cg_J_inverse_h":d_exact}
    outputs={}
    for count in (16,64):
        rate=jacobi(matrix,diag,rhs,count)
        outputs[str(count)]={
            "unprojected_volume_rate_relative_residual_l2":float(
                np.linalg.norm(fixed+matrix@rate)/np.linalg.norm(fixed)),
            "unprojected_pressure_weight_power":float(h@rate),
            "projected":{name:score(matrix,fixed,h,area,inner,rate,direction)
                         for name,direction in directions.items()},
        }
    return {
        "time":t,"initial_sha256":digest(initial_file),
        "checkpoint_sha256":digest(checkpoint),"active_cells":len(active),
        "faces_touching_active":nfaces,"cg_direction_iterations":iterations[0],
        "cg_direction_relative_residual":float(np.linalg.norm(matrix@d_exact-h)/np.linalg.norm(h)),
        "direction_jacobi16_relative_residual":float(np.linalg.norm(matrix@d16-h)/np.linalg.norm(h)),
        "direction_jacobi64_relative_residual":float(np.linalg.norm(matrix@d64-h)/np.linalg.norm(h)),
        "rate_cases":outputs,"elapsed_seconds":time.monotonic()-started,
    }


def main():
    started=time.monotonic()
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    result={
        "cases":{label:case(label) for label in ("aligned","jitter")},
        "dt_for_linear_area_prediction":DT,
        "elapsed_seconds":time.monotonic()-started,
        "scope":"Read-only finite-weight linear geometry and energy projection at native t=7.426; not a native RK4 continuation or long-time stability result.",
    }
    with OUTPUT.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    summary={label:{count:{name:row["volume_change_disk_p99"]
                            for name,row in c["projected"].items()}
                    for count,c in data["rate_cases"].items()}
             for label,data in result["cases"].items()}
    print(json.dumps({"summary_disk_p99":summary,"elapsed_seconds":result["elapsed_seconds"]},
                     indent=2,allow_nan=False))


if __name__=="__main__":
    main()

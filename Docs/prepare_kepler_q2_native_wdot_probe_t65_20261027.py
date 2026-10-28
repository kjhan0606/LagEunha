"""Prepare fixed-ID disk mask and independent W0 Jacobi rates at t=6.5."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_newton_volume_t6_20261026 import newton_jacobian
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from laguerre_weight_rate_actual_t625_20261026 import face_volume_rate
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
DOCS = Path(__file__).resolve().parent
CASES = {
    "aligned": "KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json",
    "jitter": "KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json",
}


def one_case(label):
    receipt = json.loads((DOCS / CASES[label]).read_text())
    if (receipt["time"] != 6.5 or len(receipt["checks"]) != 17
            or not all(receipt["checks"].values())):
        raise ValueError(f"Unaccepted {label} checkpoint")
    if digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]:
        raise ValueError(f"Changed {label} checkpoint")
    p0, _ = ordered(Path(receipt["initial"]), N)
    p, _ = ordered(Path(receipt["endpoint"]), N)
    xy0 = np.column_stack((p0["x"].astype(float)-3., p0["y"].astype(float)-3.))
    r0 = np.linalg.norm(xy0, axis=1)
    active = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    flag = np.zeros(N, dtype=np.uint8)
    flag[active] = 1
    lookup = np.full(N, -1, dtype=int)
    lookup[active] = np.arange(len(active))
    xy = np.column_stack((p["x"].astype(float)-3., p["y"].astype(float)-3.))
    velocity = np.column_stack((p["vx"].astype(float), p["vy"].astype(float)))
    weight = np.zeros(N)
    area, meta, triangles, centers = power_areas(
        xy, weight, active, lookup, return_geometry=True)
    if meta["hidden_active"] or np.any(~np.isfinite(area) | (area <= 0)):
        raise ValueError(f"Invalid {label} W0 geometry")
    jac, diagonal = newton_jacobian(xy, lookup, len(active), triangles, centers)
    fixed_rate, _, faces = face_volume_rate(
        xy, velocity, weight, np.zeros(N), active, lookup, triangles, centers)
    rhs = -fixed_rate
    local = np.zeros(len(active))
    for _ in range(16):
        local += .8*(rhs-jac@local)/diagonal
    rate = np.zeros(N)
    rate[active] = local
    return flag, rate, {"endpoint_sha256": receipt["endpoint_sha256"],
                        "active_cells": len(active), "faces_touching_active": faces,
                        "reference_rate_abs_max": float(np.max(np.abs(local))),
                        "reference_residual_over_rhs": float(
                            np.linalg.norm(rhs-jac@local)/np.linalg.norm(rhs))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mask", type=Path, required=True)
    ap.add_argument("--rates", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    args = ap.parse_args()
    cases = {label: one_case(label) for label in CASES}
    flag = cases["aligned"][0]
    if not np.array_equal(flag, cases["jitter"][0]):
        raise ValueError("Aligned and jittered active IDs differ")
    with args.mask.open("xb") as stream:
        flag.tofile(stream)
    with args.rates.open("xb") as stream:
        np.savez_compressed(stream, aligned=cases["aligned"][1],
                            jitter=cases["jitter"][1], active=flag)
    result = {"mask": str(args.mask.resolve()), "mask_sha256": digest(args.mask),
              "mask_bytes": args.mask.stat().st_size,
              "rates": str(args.rates.resolve()), "rates_sha256": digest(args.rates),
              "cases": {label: cases[label][2] for label in CASES},
              "scope": "Independent W0 CPU power-geometry Jacobi16 reference for one native GPU FaceCSR probe at each accepted t=6.5 restart; no weighted hydro evolution."}
    with args.manifest.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

"""Export accepted t=6.5 power faces and CPU Jacobi16 reference to CUDA."""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_newton_volume_t6_20261026 import newton_jacobian
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from laguerre_weight_rate_actual_t625_20261026 import face_volume_rate
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
MAGIC = 0x4B574452


def finite_power_edges(triangles, centers, active_flag):
    edges = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]],
                       triangles[:, [2, 0]]))
    tid = np.tile(np.arange(len(triangles)), 3)
    edges.sort(axis=1)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    edges, tid = edges[order], tid[order]
    starts = np.r_[0, 1+np.flatnonzero(np.any(edges[1:] != edges[:-1], axis=1))]
    counts = np.diff(np.r_[starts, len(edges)])
    touching = active_flag[edges[starts, 0]] | active_flag[edges[starts, 1]]
    if np.any(counts[touching] != 2):
        raise ValueError("Unpaired power face touching disk")
    starts = starts[touching]
    return (edges[starts, 0], edges[starts, 1],
            centers[tid[starts]], centers[tid[starts+1]])


def write_array(stream, array, dtype):
    np.asarray(array, dtype=dtype).tofile(stream)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accepted_receipt", type=Path)
    ap.add_argument("static_weights", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    args = ap.parse_args()
    receipt = json.loads(args.accepted_receipt.read_text())
    if (receipt["time"] != 6.5 or len(receipt["checks"]) != 17
            or not all(receipt["checks"].values())):
        raise ValueError("Unaccepted native checkpoint")
    if digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]:
        raise ValueError("Native checkpoint changed")
    p0, _ = ordered(Path(receipt["initial"]), N)
    p, _ = ordered(Path(receipt["endpoint"]), N)
    xy0 = np.column_stack((p0["x"].astype(float)-3., p0["y"].astype(float)-3.))
    xy = np.column_stack((p["x"].astype(float)-3., p["y"].astype(float)-3.))
    velocity = np.column_stack((p["vx"].astype(float), p["vy"].astype(float)))
    r0 = np.linalg.norm(xy0, axis=1)
    active = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    flag = np.zeros(N, dtype=np.uint8)
    flag[active] = 1
    lookup = np.full(N, -1, dtype=int)
    lookup[active] = np.arange(len(active))
    with np.load(args.static_weights) as saved:
        if not np.array_equal(saved["active"], active):
            raise ValueError("Static weight active set changed")
        weight = np.zeros(N)
        weight[active] = saved["weight"]
    area, metadata, triangles, centers = power_areas(
        xy, weight, active, lookup, return_geometry=True)
    if metadata["hidden_active"] or np.any(~np.isfinite(area) | (area <= 0)):
        raise ValueError("Invalid power geometry")
    jac, diagonal = newton_jacobian(xy, lookup, len(active), triangles, centers)
    fixed_rate, _, nfaces = face_volume_rate(
        xy, velocity, weight, np.zeros(N), active, lookup, triangles, centers)
    rhs = -fixed_rate
    rate = np.zeros(len(active))
    for _ in range(16):
        rate += .8*(rhs-jac@rate)/diagonal
    expected = np.zeros(N)
    expected[active] = rate
    i, j, c1, c2 = finite_power_edges(triangles, centers, flag.astype(bool))
    if len(i) != nfaces:
        raise ValueError(f"Face graph differs from Jacobian: {len(i)} != {nfaces}")
    owner = np.r_[i, j]
    neighbor = np.r_[j, i]
    corner1 = np.vstack((c1-xy[i], c2-xy[j]))
    corner2 = np.vstack((c2-xy[i], c1-xy[j]))
    order = np.lexsort((neighbor, owner))
    owner, neighbor = owner[order], neighbor[order]
    corner1, corner2 = corner1[order], corner2[order]
    offset = np.r_[0, np.cumsum(np.bincount(owner, minlength=N))]
    if offset[-1] != 2*nfaces or len(owner) != 2*nfaces:
        raise ValueError("Directed CSR face count mismatch")
    with args.output.open("xb") as stream:
        stream.write(struct.pack("<4i", MAGIC, N, len(owner), 16))
        write_array(stream, offset, "<i4")
        write_array(stream, neighbor, "<i4")
        write_array(stream, np.zeros(len(owner)), "<i4")
        for array in (corner1[:, 0], corner1[:, 1], corner2[:, 0], corner2[:, 1],
                      p["x"], p["y"], p["vx"], p["vy"]):
            write_array(stream, array, "<f8")
        write_array(stream, flag, "u1")
        write_array(stream, expected, "<f8")
    sha = digest(args.output)
    result = {"binary": str(args.output.resolve()), "sha256": sha,
              "accepted_endpoint_sha256": receipt["endpoint_sha256"],
              "static_weights_sha256": digest(args.static_weights),
              "sites": N, "active": len(active),
              "undirected_faces_touching_active": nfaces,
              "directed_face_rows": len(owner),
              "reference": "CPU Jacobi16, omega 0.8, power-area Jacobian and actual native t=6.5 velocities",
              "reference_rate_abs_max": float(np.max(np.abs(rate))),
              "scope": "Read-only exported geometry for independent CUDA operator check; no native RK4 evolution."}
    with args.manifest.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("sha256", "active", "directed_face_rows", "reference_rate_abs_max")}, indent=2))


if __name__ == "__main__":
    main()

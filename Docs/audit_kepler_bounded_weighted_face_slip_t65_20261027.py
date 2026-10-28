"""Compare geometric and Kepler-material normal speeds on disk-disk faces."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
EPS = 0.01175


def finite_faces(triangles, centers, active_flag):
    edges = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]],
                       triangles[:, [2, 0]]))
    triangle_ids = np.tile(np.arange(len(triangles)), 3)
    edges.sort(axis=1)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    edges, triangle_ids = edges[order], triangle_ids[order]
    starts = np.r_[0, 1+np.flatnonzero(np.any(edges[1:] != edges[:-1], axis=1))]
    counts = np.diff(np.r_[starts, len(edges)])
    touching = active_flag[edges[starts, 0]] | active_flag[edges[starts, 1]]
    if np.any(counts[touching] != 2):
        raise ValueError("Non-finite or nonmanifold face touching the disk")
    starts = starts[(counts == 2) & active_flag[edges[starts, 0]]
                    & active_flag[edges[starts, 1]]]
    i, j = edges[starts, 0], edges[starts, 1]
    c1, c2 = centers[triangle_ids[starts]], centers[triangle_ids[starts+1]]
    midpoint = .5*(c1+c2)
    length = np.linalg.norm(c2-c1, axis=1)
    if not np.all(np.isfinite(length) & (length > 0)):
        raise ValueError("Invalid finite face length")
    return i, j, midpoint, length


def metrics(length, slip):
    abs_slip = np.abs(slip)
    return {"faces": len(slip), "total_face_length": float(np.sum(length, dtype=np.longdouble)),
            "length_weighted_mean_absolute_slip": float(
                np.sum(length*abs_slip, dtype=np.longdouble)
                / np.sum(length, dtype=np.longdouble)),
            "length_weighted_rms_slip": float(np.sqrt(
                np.sum(length*slip*slip, dtype=np.longdouble)
                / np.sum(length, dtype=np.longdouble))),
            "integrated_absolute_slip": float(np.sum(length*abs_slip, dtype=np.longdouble)),
            "maximum_absolute_slip": float(np.max(abs_slip))}


def diagram(xy, velocity, r0, active, lookup, weight, rate):
    _, meta, triangles, centers = power_areas(
        xy, weight, active, lookup, return_geometry=True)
    if meta["hidden_active"]:
        raise ValueError("Hidden active power cell")
    active_flag = lookup >= 0
    i, j, c, length = finite_faces(triangles, centers, active_flag)
    dr = xy[j]-xy[i]
    d = np.linalg.norm(dr, axis=1)
    normal = dr/d[:, None]
    fixed = (np.einsum("ij,ij->i", xy[j], velocity[j])
             - np.einsum("ij,ij->i", xy[i], velocity[i])
             - np.einsum("ij,ij->i", velocity[j]-velocity[i], c))/d
    moving = fixed+(rate[i]-rate[j])/(2*d)
    rface2 = np.einsum("ij,ij->i", c, c)
    material_velocity = np.column_stack((-c[:, 1], c[:, 0]))
    material_velocity *= (rface2+EPS*EPS)[:, None]**(-.75)
    material = np.einsum("ij,ij->i", material_velocity, normal)
    same_shell = np.abs(r0[i]-r0[j]) < 1e-10
    return {"geometry": meta,
            "all_disk_disk": {"fixed": metrics(length, fixed-material),
                              "moving": metrics(length, moving-material)},
            "same_initial_shell": {"fixed": metrics(length[same_shell], (fixed-material)[same_shell]),
                                   "moving": metrics(length[same_shell], (moving-material)[same_shell])},
            "different_initial_shell": {"fixed": metrics(length[~same_shell], (fixed-material)[~same_shell]),
                                        "moving": metrics(length[~same_shell], (moving-material)[~same_shell])}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accepted_receipt", type=Path)
    ap.add_argument("static_weights", type=Path)
    ap.add_argument("weight_rates", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    receipt = json.loads(args.accepted_receipt.read_text())
    if (receipt["time"] != 6.5 or len(receipt["checks"]) != 17
            or not all(receipt["checks"].values())):
        raise ValueError("Unaccepted t=6.5 snapshot")
    if (digest(Path(receipt["initial"])) != receipt["initial_sha256"]
            or digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]):
        raise ValueError("Accepted snapshots changed")
    p0, _ = ordered(Path(receipt["initial"]), N)
    p, _ = ordered(Path(receipt["endpoint"]), N)
    xy0 = np.column_stack((p0["x"].astype(float)-3., p0["y"].astype(float)-3.))
    xy = np.column_stack((p["x"].astype(float)-3., p["y"].astype(float)-3.))
    velocity = np.column_stack((p["vx"].astype(float), p["vy"].astype(float)))
    r0 = np.linalg.norm(xy0, axis=1)
    active = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    lookup = np.full(N, -1, dtype=int)
    lookup[active] = np.arange(len(active))
    with np.load(args.static_weights) as saved:
        if not np.array_equal(saved["active"], active):
            raise ValueError("Static active set differs")
        local_weight = saved["weight"].copy()
    with np.load(args.weight_rates) as saved:
        if not np.array_equal(saved["active"], active):
            raise ValueError("Rate active set differs")
        local_rate = saved["weight_rate"].copy()
        local_rate16 = saved["rate_16"].copy()
    weight = np.zeros(N)
    weight[active] = local_weight
    rate = np.zeros(N)
    rate[active] = local_rate
    rate16 = np.zeros(N)
    rate16[active] = local_rate16
    zero = np.zeros(N)
    voronoi = diagram(xy, velocity, r0, active, lookup, zero, zero)
    weighted_fixed = diagram(xy, velocity, r0, active, lookup, weight, zero)
    weighted_full = diagram(xy, velocity, r0, active, lookup, weight, rate)
    weighted_16 = diagram(xy, velocity, r0, active, lookup, weight, rate16)
    result = {"accepted_receipt": str(args.accepted_receipt.resolve()),
              "endpoint_sha256": receipt["endpoint_sha256"],
              "static_weights_sha256": digest(args.static_weights),
              "weight_rates_sha256": digest(args.weight_rates),
              "time": 6.5, "active_disk_sites": len(active),
              "voronoi_W0": voronoi,
              "static_power_Wdot0": weighted_fixed,
              "static_power_Wdot16": weighted_16,
              "static_power_Wdot_full": weighted_full,
              "scope": "Analytic softened Kepler material speed at each disk-disk face centroid versus geometric face speed. A volume-rate solve constrains the signed sum, not every face. This is not native mass flux or evolved weighted hydro."}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    for label in ("voronoi_W0", "static_power_Wdot0", "static_power_Wdot16", "static_power_Wdot_full"):
        row = result[label]["all_disk_disk"]["moving"]
        print(label, row["faces"], row["length_weighted_mean_absolute_slip"])


if __name__ == "__main__":
    main()

"""One-time IC quadrature: compare generator and cell-mean central gravity.

This is an offline spatial consistency diagnostic, not a replacement for the
native point-source RK4 evolution or a proposed per-step gravity integral.
"""

import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import Voronoi, cKDTree


BASE = Path("/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026")
DOCS = Path(__file__).resolve().parent
EPS = .01175
R0 = .10080469
GLASSES = {
    "full_q2": (BASE / "q2/glass.bin",
                "3ea6edbf56f882d648abea0bef620dac1564f9146194f332f2566af4412ad4a4", None),
    "narrow": (BASE / "q2_profile_buffered_t001/glass.bin",
               "03448f6151f85195b49e8ea7bd7543611094b8b0a957c3891134d1d83f0156a7",
               DOCS / "KEPLER_Q2_PROFILE_BUFFERED_PERIODIC_AREAS_2026-10-27.bin"),
    "wide": (BASE / "q2_profile_wide_t001/glass.bin",
             "de50b97639a243f9384fc7b4bfc402f9eccb52386f30ecdd3450f7b0acaa9e8f",
             DOCS / "KEPLER_Q2_PROFILE_WIDE_PERIODIC_AREAS_2026-10-27.bin"),
}


def load_glass(path, expected):
    raw = path.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    if sha != expected:
        raise ValueError(f"Unexpected glass SHA: {path}")
    n, lx, ly = struct.unpack_from("<idd", raw)
    if (lx, ly) != (6., 6.) or len(raw) != 20 + 16*n:
        raise ValueError("Invalid native glass")
    x = np.frombuffer(raw, "<f8", count=n, offset=20)
    y = np.frombuffer(raw, "<f8", count=n, offset=20 + 8*n)
    return np.column_stack((x, y)), sha


def load_areas(path, n):
    raw = path.read_bytes()
    count, lx, ly = struct.unpack_from("<idd", raw)
    if (count, lx, ly) != (n, 6., 6.) or len(raw) != 20 + 8*n:
        raise ValueError("Invalid independent periodic area sidecar")
    return np.frombuffer(raw, "<f8", count=n, offset=20), hashlib.sha256(raw).hexdigest()


def local_cells(xy, chosen):
    radius = np.linalg.norm(xy - 3., axis=1)
    local_id = np.flatnonzero(radius < .4)
    outside_id = np.flatnonzero(radius >= .4)
    vor = Voronoi(xy[local_id], qhull_options="Qbb Qc Qz")
    excluded = cKDTree(xy[outside_id])
    reverse = {int(global_id): k for k, global_id in enumerate(local_id)}
    result = {}
    for global_id in chosen:
        site = xy[global_id]
        region = vor.regions[vor.point_region[reverse[int(global_id)]]]
        if not region or -1 in region:
            raise ValueError("Unbounded central cell")
        vertices = vor.vertices[region]
        nearest_excluded, _ = excluded.query(vertices)
        own_distance = np.linalg.norm(vertices-site, axis=1)
        if not np.all(nearest_excluded > own_distance + 1e-10):
            raise ValueError("The local Voronoi crop missed a relevant site")
        result[int(global_id)] = vertices
    return result


def integrate_cell(site, vertices, order):
    q, w = np.polynomial.legendre.leggauss(order)
    q, w = (q+1.)/2., w/2.
    u, v = q[:, None, None], q[None, :, None]
    weights = w[:, None] * w[None, :] * q[:, None]
    area = mass = 0.
    force = np.zeros(2)
    weighted_force = np.zeros(2)
    for a, b in zip(vertices, np.roll(vertices, -1, axis=0)):
        det = abs(np.linalg.det(np.stack((a-site, b-site))))
        point = site + u*(a-site) + u*v*(b-a)
        relative = point - 3.
        r2 = np.sum(relative*relative, axis=2)
        gravity = -relative / (r2[:, :, None] + EPS*EPS)**1.5
        rho = .01 + (np.sqrt(r2)/.5)**3
        dA = weights*det
        area += float(np.sum(dA))
        mass += float(np.sum(rho*dA))
        force += np.sum(gravity*dA[:, :, None], axis=(0, 1))
        weighted_force += np.sum(gravity*rho[:, :, None]*dA[:, :, None], axis=(0, 1))
    g_site = -(site-3.) / (np.dot(site-3., site-3.) + EPS*EPS)**1.5
    radial = g_site / np.linalg.norm(g_site)
    rho_site = .01 + (np.linalg.norm(site-3.)/.5)**3
    return {
        "area": area,
        "area_mean_force_radial_over_generator": float(np.dot(force/area, radial)/np.linalg.norm(g_site)),
        "profile_weighted_force_radial_over_generator": float(np.dot(weighted_force/mass, radial)/np.linalg.norm(g_site)),
        "profile_cell_average_density_over_generator_density": float(mass/area/rho_site),
        "area_mean_force_tangential_over_generator": float(np.linalg.det(np.stack((force/area, radial)))/np.linalg.norm(g_site)),
    }


def summarize(rows):
    keys = ("area_mean_force_radial_over_generator",
            "profile_weighted_force_radial_over_generator",
            "profile_cell_average_density_over_generator_density")
    return {key: {"min": float(min(row[key] for row in rows)),
                  "mean": float(np.mean([row[key] for row in rows])),
                  "max": float(max(row[key] for row in rows))} for key in keys}


def main():
    output = DOCS / "KEPLER_SPARSE_CORE_GRAVITY_IC_2026-10-27.json"
    if output.exists():
        raise FileExistsError(output)
    datasets = {key: load_glass(path, sha) for key, (path, sha, _) in GLASSES.items()}
    full_xy = datasets["full_q2"][0]
    full_tree = cKDTree(full_xy)
    report = {"softening": EPS, "local_voronoi_crop_radius": .4,
              "quadrature_orders": [32, 64], "cases": {}}
    for name in ("narrow", "wide"):
        xy, glass_sha = datasets[name]
        r = np.linalg.norm(xy-3., axis=1)
        first = np.flatnonzero(abs(r-R0) < 1e-6)
        if len(first) != 8:
            raise ValueError(f"Expected eight first-ring cells, got {len(first)}")
        polygon = local_cells(xy, first)
        sidecar, area_sha = load_areas(GLASSES[name][2], len(xy))
        rows = []
        for idx in first:
            a = integrate_cell(xy[idx], polygon[int(idx)], 32)
            b = integrate_cell(xy[idx], polygon[int(idx)], 64)
            b.update({"site_id": int(idx), "radius": float(r[idx]),
                      "sidecar_area": float(sidecar[idx]),
                      "area_relative_to_sidecar": float(b["area"]/sidecar[idx]-1.),
                      "quadrature_force_ratio_difference_32_to_64": float(
                          b["profile_weighted_force_radial_over_generator"]-
                          a["profile_weighted_force_radial_over_generator"])})
            if abs(b["area"]/sidecar[idx]-1.) > 1e-8:
                raise ValueError("Local central polygon disagrees with periodic sidecar")
            rows.append(b)
        report["cases"][name] = {"glass_sha256": glass_sha, "area_sha256": area_sha,
                                  "first_ring_cells": len(rows), "summary": summarize(rows),
                                  "maximum_area_relative_error": float(max(abs(row["area_relative_to_sidecar"]) for row in rows)),
                                  "maximum_32_to_64_force_ratio_difference": float(max(abs(row["quadrature_force_ratio_difference_32_to_64"]) for row in rows)),
                                  "cells": rows}
    sparse_xy = datasets["narrow"][0]
    sparse_r = np.linalg.norm(sparse_xy-3., axis=1)
    first = np.flatnonzero(abs(sparse_r-R0) < 1e-6)
    _, full_nearest = full_tree.query(sparse_xy[first])
    if len(set(full_nearest)) != 8:
        raise ValueError("Full-q2 comparison selected duplicate cells")
    full_cells = local_cells(full_xy, full_nearest)
    rows = []
    for idx in full_nearest:
        row = integrate_cell(full_xy[idx], full_cells[int(idx)], 64)
        row.update({"site_id": int(idx), "radius": float(np.linalg.norm(full_xy[idx]-3.))})
        rows.append(row)
    report["cases"]["full_q2_near_first_ring"] = {
        "glass_sha256": datasets["full_q2"][1], "sampled_cells": len(rows),
        "selection": "Nearest full-q2 site to each of the narrow sparse first-ring sites",
        "summary": summarize(rows), "cells": rows}
    report["interpretation_limit"] = (
        "Measures initial spatial point-force versus cell-average mismatch. "
        "It does not prove native heating, and does not prescribe per-step cell-integrated gravity.")
    with output.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps({key: {"summary": value["summary"],
                            "count": value.get("first_ring_cells", value.get("sampled_cells"))}
                      for key, value in report["cases"].items()}, indent=2))


if __name__ == "__main__":
    main()

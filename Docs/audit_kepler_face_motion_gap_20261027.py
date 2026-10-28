"""Measure exact-orbit material versus Voronoi face-normal motion.

This is a geometric kinematic probe, not a native HLLC flux or heat audit.
Only faces incident on initial disk cells are used; these are bounded well
inside the 6x6 box and do not depend on periodic images at t=0.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi

from audit_kepler_polar_checkpoint_20261026 import CX, CY, EPS, load_ordered


def weighted_mean(value, weight):
    return float(np.sum(value * weight, dtype=np.longdouble) /
                 np.sum(weight, dtype=np.longdouble))


def orbital_velocity(xy):
    centered = xy - (CX, CY)
    r2 = np.sum(centered * centered, axis=1)
    omega = (r2 + EPS * EPS) ** (-.75)
    return np.column_stack((-centered[:, 1] * omega,
                             centered[:, 0] * omega))


def ballistic_xy(xy, time):
    centered = xy - (CX, CY)
    r2 = np.sum(centered * centered, axis=1)
    angle = time * (r2 + EPS * EPS) ** (-.75)
    ca, sa = np.cos(angle), np.sin(angle)
    return np.column_stack((CX + ca * centered[:, 0] - sa * centered[:, 1],
                            CY + sa * centered[:, 0] + ca * centered[:, 1]))


def selected_areas(vor, selected):
    area = np.empty(len(selected), dtype=float)
    for k, i in enumerate(selected):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded disk polygon: {i}")
        polygon = vor.vertices[region]
        area[k] = .5 * abs(np.dot(polygon[:, 0], np.roll(polygon[:, 1], -1)) -
                             np.dot(polygon[:, 1], np.roll(polygon[:, 0], -1)))
    if not np.all(np.isfinite(area) & (area > 0)):
        raise ValueError("Nonpositive or nonfinite disk area")
    return area


def face_summary(mask, length, gap, geom, material):
    if not np.any(mask):
        raise ValueError("Empty face group")
    w = length[mask]
    g = gap[mask]
    return {
        "faces": int(np.sum(mask)),
        "total_length": float(np.sum(w)),
        "mean_abs_material_minus_geometric_normal_speed": weighted_mean(abs(g), w),
        "rms_material_minus_geometric_normal_speed":
            float(np.sqrt(weighted_mean(g * g, w))),
        "p95_abs_material_minus_geometric_normal_speed":
            float(np.quantile(abs(g), .95)),
        "length_integral_abs_gap": float(np.sum(w * abs(g))),
        "length_integral_gap_squared": float(np.sum(w * g * g)),
        "mean_abs_geometric_normal_speed": weighted_mean(abs(geom[mask]), w),
        "mean_abs_material_normal_speed": weighted_mean(abs(material[mask]), w),
    }


def audit(path):
    p, meta = load_ordered(path)
    if float(meta["t"]) != 0 or not np.all(p["w2"] == 0):
        raise ValueError("Expected complete native t=0, W=0 state")
    xy = np.column_stack((p["x"], p["y"]))
    centered = xy - (CX, CY)
    r = np.linalg.norm(centered, axis=1)
    disk = (r >= .5) & (r < 2.)
    velocity = orbital_velocity(xy)
    vor = Voronoi(xy, qhull_options="Qbb Qc Qz")
    all_pairs = vor.ridge_points
    incident = disk[all_pairs[:, 0]] | disk[all_pairs[:, 1]]
    pairs = all_pairs[incident]
    a, b = pairs[:, 0], pairs[:, 1]
    ridges = [vertices for vertices, selected in
              zip(vor.ridge_vertices, incident) if selected]
    if any(len(vertices) != 2 or min(vertices) < 0 for vertices in ridges):
        raise ValueError("An initial disk face is unbounded")
    ridge_vertices = np.asarray(ridges, dtype=np.intp)
    corners = vor.vertices[ridge_vertices]
    length = np.linalg.norm(corners[:, 1] - corners[:, 0], axis=1)
    distance_vec = xy[b] - xy[a]
    distance = np.linalg.norm(distance_vec, axis=1)
    if not np.all(np.isfinite(length) & (length > 0)) or not np.all(distance > 0):
        raise ValueError("Degenerate initial disk face")
    normal = distance_vec / distance[:, None]
    centroid = .5 * (corners[:, 0] + corners[:, 1])
    mid_generators = .5 * (xy[a] + xy[b])
    mean_velocity = .5 * (velocity[a] + velocity[b])
    delta_velocity = velocity[b] - velocity[a]
    # Exact W=0 bisector normal speed at the face centroid. It matches the
    # native voronoi_face_rotation expression with unsteered generators.
    geom_speed = (np.einsum("ij,ij->i", mean_velocity, normal)
                  - np.einsum("ij,ij->i", delta_velocity,
                              centroid - mid_generators) / distance)
    material_speed = np.einsum("ij,ij->i", orbital_velocity(centroid), normal)
    gap = material_speed - geom_speed
    if not np.all(np.isfinite(gap)):
        raise ValueError("Nonfinite face-motion gap")

    both_disk = disk[a] & disk[b]
    same_shell = both_disk & (np.round(r[a], 10) == np.round(r[b], 10))
    cross_shell = both_disk & ~same_shell
    disk_guard = ~both_disk
    groups = {
        "all_disk_incident": face_summary(np.ones(len(a), dtype=bool), length,
                                          gap, geom_speed, material_speed),
        "disk_disk_same_shell": face_summary(same_shell, length, gap,
                                             geom_speed, material_speed),
        "disk_disk_cross_shell": face_summary(cross_shell, length, gap,
                                              geom_speed, material_speed),
        "disk_guard": face_summary(disk_guard, length, gap,
                                   geom_speed, material_speed),
    }

    # Signed face fluxes use the outward normal from a to b. Their sum over
    # each complete disk polygon gives its geometric and material area rates.
    geom_rate = np.zeros(len(p), dtype=float)
    material_rate = np.zeros(len(p), dtype=float)
    np.add.at(geom_rate, a, length * geom_speed)
    np.add.at(geom_rate, b, -length * geom_speed)
    np.add.at(material_rate, a, length * material_speed)
    np.add.at(material_rate, b, -length * material_speed)
    selected = np.flatnonzero(disk)
    area = selected_areas(vor, selected)
    # Independent finite-difference geometry check. These Voronoi diagrams
    # are reconstructed anew; their cell areas do not use the face rates.
    h = 1e-4
    area_plus = selected_areas(Voronoi(ballistic_xy(xy, h),
                                       qhull_options="Qbb Qc Qz"), selected)
    area_minus = selected_areas(Voronoi(ballistic_xy(xy, -h),
                                        qhull_options="Qbb Qc Qz"), selected)
    finite_difference_rate = (area_plus - area_minus) / (2 * h)
    rate_error = geom_rate[selected] - finite_difference_rate
    mass = p["mass"][selected].astype(float)
    geom_relative = geom_rate[selected] / area
    material_relative = material_rate[selected] / area
    difference_relative = material_relative - geom_relative
    result = {
        "snapshot": str(path),
        "disk_cells": int(np.sum(disk)),
        "groups": groups,
        "disk_cell_relative_area_rate": {
            "geometric_abs_mass_weighted": weighted_mean(abs(geom_relative), mass),
            "material_midpoint_abs_mass_weighted": weighted_mean(abs(material_relative), mass),
            "material_minus_geometric_abs_mass_weighted":
                weighted_mean(abs(difference_relative), mass),
            "material_minus_geometric_p95_abs":
                float(np.quantile(abs(difference_relative), .95)),
            "geometric_signed_area_sum": float(np.sum(geom_rate[selected])),
        },
        "independent_geometric_rate_check": {
            "central_difference_dt": h,
            "native_initial_area_max_relative_difference":
                float(np.max(abs(area / p["volume"][selected].astype(float) - 1))),
            "mean_abs_rate_error_relative_to_cell_area":
                weighted_mean(abs(rate_error) / area, mass),
            "p95_abs_rate_error_relative_to_cell_area":
                float(np.quantile(abs(rate_error) / area, .95)),
            "max_abs_rate_error_relative_to_cell_area":
                float(np.max(abs(rate_error) / area)),
        },
        "scope": "Analytic Kepler material speed minus exact W0 geometric face speed at native t=0 polygon face centroids. Not native HLLC contact speed, mass flux, entropy, or heating.",
    }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("snapshots", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = [audit(path) for path in args.snapshots]
    if len(result) == 2:
        a, b = (row["groups"]["disk_disk_cross_shell"] for row in result)
        result.append({"second_over_first_cross_shell": {
            key: b[key] / a[key] for key in
            ("mean_abs_material_minus_geometric_normal_speed",
             "rms_material_minus_geometric_normal_speed",
             "length_integral_abs_gap", "length_integral_gap_squared")}})
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

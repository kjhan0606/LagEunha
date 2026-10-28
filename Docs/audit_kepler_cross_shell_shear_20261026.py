"""Measure geometric shear exposure on native W=0 Kepler disk faces.

This is a read-only kinematic diagnostic, not a viscous-heating measurement.
It uses the unweighted Voronoi faces of complete native t=0 snapshots.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi

from audit_kepler_polar_checkpoint_20261026 import CX, CY, EPS, load_ordered


def group_stats(mask, length, separation, shear, approach, radial_normal):
    weights = length[mask]
    values = {
        "faces": int(np.sum(mask)),
        "face_length": float(np.sum(weights)),
        "length_weighted_center_distance": float(np.average(separation[mask], weights=weights)),
        "length_weighted_local_shear_speed": float(np.average(shear[mask], weights=weights)),
        "length_integral_local_shear_speed": float(np.sum(weights * shear[mask])),
        "length_integral_local_shear_speed_squared": float(np.sum(weights * shear[mask] ** 2)),
        "length_weighted_approach_speed": float(np.average(approach[mask], weights=weights)),
        "length_integral_approach_speed": float(np.sum(weights * approach[mask])),
        "length_weighted_radial_normal_abs": float(np.average(radial_normal[mask], weights=weights)),
    }
    return values


def audit(path):
    p, meta = load_ordered(path)
    if float(meta["t"]) != 0 or not np.all(p["w2"] == 0):
        raise ValueError("Expected complete native t=0 W=0 snapshot")
    xy = np.column_stack((p["x"], p["y"]))
    centered = xy - np.array((CX, CY))
    # Match the established native face audit's disk-boundary convention.
    radius = np.hypot(centered[:, 0], centered[:, 1])
    disk = (radius >= .5) & (radius < 2.)
    vor = Voronoi(xy, qhull_options="Qbb Qc Qz")
    all_pairs = vor.ridge_points
    keep = disk[all_pairs[:, 0]] & disk[all_pairs[:, 1]]
    pairs = all_pairs[keep]
    a, b = pairs[:, 0], pairs[:, 1]
    ridge_vertices = np.asarray(
        [vertices for vertices, selected in zip(vor.ridge_vertices, keep) if selected],
        dtype=np.intp,
    )
    if ridge_vertices.ndim != 2 or ridge_vertices.shape[1] != 2 or np.any(ridge_vertices < 0):
        raise ValueError("Disk contains an unbounded Voronoi face")
    length = np.linalg.norm(
        vor.vertices[ridge_vertices[:, 0]] - vor.vertices[ridge_vertices[:, 1]], axis=1
    )
    separation_vec = xy[b] - xy[a]
    separation = np.linalg.norm(separation_vec, axis=1)
    if np.any(length <= 0) or np.any(separation <= 0):
        raise ValueError("Degenerate disk face")
    face_normal = separation_vec / separation[:, None]

    # In the co-rotating local frame, r_mid * Delta Omega is the relative
    # azimuthal drift of adjacent orbital rings across this face.
    omega = (radius * radius + EPS * EPS) ** (-.75)
    shear = .5 * (radius[a] + radius[b]) * np.abs(omega[a] - omega[b])
    velocity = np.column_stack((p["vx"], p["vy"]))
    pair_normal_speed = np.sum((velocity[b] - velocity[a]) * face_normal, axis=1)
    approach = np.maximum(-pair_normal_speed, 0.)
    midpoint = .5 * (centered[a] + centered[b])
    midpoint_radius = np.linalg.norm(midpoint, axis=1)
    radial_normal = np.abs(np.sum(face_normal * midpoint, axis=1)) / midpoint_radius

    shell_radius, shell = np.unique(np.round(radius, 10), return_inverse=True)
    same_shell = shell[a] == shell[b]
    cross_shell = ~same_shell
    if not np.any(cross_shell) or not np.any(same_shell):
        raise ValueError("Expected same-shell and cross-shell faces")
    return {
        "snapshot": str(path),
        "time": float(meta["t"]),
        "disk_cells": int(np.sum(disk)),
        "polar_ring_radii": int(len(np.unique(shell[disk]))),
        "same_shell": group_stats(same_shell, length, separation, shear, approach, radial_normal),
        "cross_shell": group_stats(cross_shell, length, separation, shear, approach, radial_normal),
        "scope": "Exact t=0 W=0 face geometry and native velocities. Local shear is r_mid*|Omega_i-Omega_j|. Approach is max[-(v_j-v_i) dot n_ij,0]. Neither quantity is measured heating or a face flux.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("snapshots", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = [audit(path) for path in args.snapshots]
    if len(result) == 2:
        a, b = (record["cross_shell"] for record in result)
        result.append({
            "comparison_second_over_first": {
                key: b[key] / a[key]
                for key in ("face_length", "length_weighted_local_shear_speed",
                            "length_integral_local_shear_speed",
                            "length_integral_local_shear_speed_squared",
                            "length_integral_approach_speed")
            }
        })
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

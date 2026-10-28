"""Read-only native Kepler polar face and shell-separation audit.

Compare complete t=0 and early native snapshots with an independent exact-orbit
Voronoi mesh. The disk is interior to the 6x6 box, so periodic copies cannot
affect its finite face segments at these times.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi

from audit_kepler_polar_checkpoint_20261026 import (
    CX, CY, EPS, GAMMA_MINUS_ONE, N, load_ordered, weighted_mean,
)


def face_graph(xy, disk, shell, same_chord):
    vor = Voronoi(xy, qhull_options="Qbb Qc Qz")
    pair = vor.ridge_points
    keep = disk[pair[:, 0]] & disk[pair[:, 1]]
    selected = pair[keep]
    a, b = selected[:, 0], selected[:, 1]
    ridges = [v for v, accept in zip(vor.ridge_vertices, keep) if accept]
    if any(len(v) != 2 or min(v) < 0 for v in ridges):
        raise ValueError("Unbounded or non-segment disk face")
    vertices = np.asarray(ridges, dtype=np.intp)
    length = np.linalg.norm(vor.vertices[vertices[:, 0]] -
                            vor.vertices[vertices[:, 1]], axis=1)
    if not np.all(np.isfinite(length) & (length > 0)):
        raise ValueError("Invalid disk face length")
    center_distance = np.linalg.norm(xy[a] - xy[b], axis=1)
    # Use the larger chord: ratio>=1 proves the pair is farther than either
    # endpoint's nearest same-shell generator.
    ratio = center_distance / np.maximum(same_chord[a], same_chord[b])
    same = shell[a] == shell[b]
    adjacent = np.abs(shell[a] - shell[b]) == 1
    cross = ~same
    keys = np.minimum(a, b).astype(np.uint64) * np.uint64(N) + np.maximum(a, b)
    return {
        "keys": np.sort(keys),
        "faces": int(len(a)),
        "same_shell_faces": int(np.sum(same)),
        "adjacent_shell_faces": int(np.sum(adjacent)),
        "other_cross_shell_faces": int(np.sum(cross & ~adjacent)),
        "face_length_sum": float(np.sum(length)),
        "same_shell_face_length_sum": float(np.sum(length[same])),
        "cross_shell_face_length_sum": float(np.sum(length[cross])),
        "cross_shell_face_length_fraction": float(np.sum(length[cross]) / np.sum(length)),
        "adjacent_shell_face_length_fraction": float(np.sum(length[adjacent]) / np.sum(length)),
        "cross_shell_center_closer_than_same_shell_chord": int(np.sum(ratio[cross] < 1)),
        "cross_shell_center_distance_over_chord_min": float(np.min(ratio[cross])),
        "cross_shell_center_distance_over_chord_p01": float(np.quantile(ratio[cross], .01)),
        "cross_shell_center_distance_over_chord_median": float(np.median(ratio[cross])),
    }


def overlap(a, b):
    common = len(np.intersect1d(a["keys"], b["keys"], assume_unique=True))
    return {
        "common_faces": common,
        "lost_from_first": len(a["keys"]) - common,
        "new_in_second": len(b["keys"]) - common,
        "jaccard": common / (len(a["keys"]) + len(b["keys"]) - common),
    }


def report(initial, endpoint):
    p0, meta0 = load_ordered(initial)
    pt, metat = load_ordered(endpoint)
    t = float(metat["t"])
    if float(meta0["t"]) != 0 or not 0 < t <= .825 + 1e-12:
        raise ValueError("Expected complete t=0 and early endpoint snapshots")
    if not np.array_equal(p0["mass"], pt["mass"]) or not np.all(pt["w2"] == 0):
        raise ValueError("Expected fixed-mass, W=0 snapshots")
    r0 = np.hypot(p0["x"] - CX, p0["y"] - CY)
    disk = (r0 >= .5) & (r0 < 2.)
    polar = r0 < 2.35
    radii, shell_local, counts = np.unique(
        np.round(r0[polar], 10), return_inverse=True, return_counts=True)
    shell = np.full(N, -1, dtype=np.intp)
    shell[polar] = shell_local
    if not np.all(shell[disk] >= 0):
        raise ValueError("Disk cell absent from polar rings")
    same_chord = np.ones(N, dtype=np.float64)
    same_chord[polar] = 2 * radii[shell_local] * np.sin(np.pi / counts[shell_local])
    if not np.all(same_chord[disk] > 0):
        raise ValueError("Invalid same-shell chord")

    xy0 = np.column_stack((p0["x"], p0["y"]))
    x = p0["x"] - CX
    y = p0["y"] - CY
    angle = t * (r0*r0 + EPS*EPS)**(-.75)
    ballistic = np.column_stack((CX + x*np.cos(angle) - y*np.sin(angle),
                                  CY + x*np.sin(angle) + y*np.cos(angle)))
    native = np.column_stack((pt["x"], pt["y"]))
    g0 = face_graph(xy0, disk, shell, same_chord)
    gb = face_graph(ballistic, disk, shell, same_chord)
    gn = face_graph(native, disk, shell, same_chord)

    entropy_log = (np.log(pt["ie"].astype(float) / p0["ie"].astype(float)) +
                   GAMMA_MINUS_ONE *
                   np.log(pt["volume"].astype(float) / p0["volume"].astype(float)))
    if not np.all(np.isfinite(entropy_log[disk])):
        raise ValueError("Nonfinite disk entropy change")
    weights = p0["mass"].astype(float)
    disk_entropy = entropy_log[disk]
    result = {
        "initial_snapshot": str(initial),
        "endpoint": str(endpoint),
        "time": t,
        "disk_cells": int(np.sum(disk)),
        "polar_rings": int(len(radii)),
        "initial": {k: v for k, v in g0.items() if k != "keys"},
        "ballistic": {k: v for k, v in gb.items() if k != "keys"},
        "native": {k: v for k, v in gn.items() if k != "keys"},
        "initial_to_ballistic_graph": overlap(g0, gb),
        "initial_to_native_graph": overlap(g0, gn),
        "ballistic_to_native_graph": overlap(gb, gn),
        "disk_entropy_log_mass_weighted_signed": weighted_mean(disk_entropy, weights[disk]),
        "disk_entropy_log_mass_weighted_abs": weighted_mean(abs(disk_entropy), weights[disk]),
        "disk_entropy_log_p95_abs": float(np.quantile(abs(disk_entropy), .95)),
        "scope": "Unweighted Voronoi face geometry and native cell entropy; no per-face flux attribution",
    }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("endpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = report(args.initial, args.endpoint)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

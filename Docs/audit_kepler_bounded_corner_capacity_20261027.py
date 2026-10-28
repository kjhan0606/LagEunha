"""Independent bounded Voronoi check of a native accepted Kepler checkpoint."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi

from audit_kepler_sparse_area_t0825_20261027 import ordered
from make_kepler_bounded_polygon_ic_20261027 import (
    BOUNDARY_RADIUS, BOUNDARY_SIDES, CENTER, GUARD_SIDES,
    clip_halfplane, polygon_area,
)
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N = 55351
DOMAIN_AREA = 0.5 * BOUNDARY_SIDES * BOUNDARY_RADIUS**2 * np.sin(
    2 * np.pi / BOUNDARY_SIDES)


def reconstruct(xy, guard_radius):
    guard_theta = np.arange(GUARD_SIDES) * 2 * np.pi / GUARD_SIDES
    guard = CENTER + guard_radius * np.column_stack(
        (np.cos(guard_theta), np.sin(guard_theta)))
    vor = Voronoi(np.vstack((xy, guard)), qhull_options="Qbb Qc Qz")
    theta = (np.arange(BOUNDARY_SIDES) + 0.5) * 2 * np.pi / BOUNDARY_SIDES
    normals = np.column_stack((np.cos(theta), np.sin(theta)))
    apothem = BOUNDARY_RADIUS * np.cos(np.pi / BOUNDARY_SIDES)
    thresholds = normals @ CENTER + apothem
    areas = np.empty(len(xy))
    corners = np.empty(len(xy), dtype=int)
    clipped = np.zeros(len(xy), dtype=bool)
    for i in range(len(xy)):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded original cell {i}")
        poly = vor.vertices[region]
        if np.max(np.linalg.norm(poly - CENTER, axis=1)) > apothem:
            clipped[i] = True
            for normal, threshold in zip(normals, thresholds):
                poly = clip_halfplane(poly, normal, threshold)
        areas[i] = polygon_area(poly)
        corners[i] = len(poly)
    if not np.all(np.isfinite(areas) & (areas > 0)):
        raise ValueError("Nonpositive reconstructed area")
    return areas, corners, clipped


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    t = snapshot_time(args.snapshot, N)
    state, _ = ordered(args.snapshot, N)
    xy = np.column_stack((state["x"].astype(float),
                          state["y"].astype(float)))
    a4, c4, b4 = reconstruct(xy, 4.0)
    a5, c5, b5 = reconstruct(xy, 5.0)
    native = state["volume"].astype(float)
    native_rel = np.abs(native / a5 - 1)
    guard_rel = np.abs(a4 / a5 - 1)
    checks = {
        "independent_guards_agree": bool(np.max(guard_rel) < 1e-8),
        "boundary_cells_agree": bool(np.array_equal(b4, b5)),
        "final_corner_counts_agree": bool(np.array_equal(c4, c5)),
        "areas_tile_polygon": bool(abs(np.sum(a5, dtype=np.longdouble)
                                       - DOMAIN_AREA) < 1e-7),
        "native_areas_match": bool(np.max(native_rel) < 1e-5),
    }
    result = {
        "snapshot": str(args.snapshot.resolve()),
        "sha256": digest(args.snapshot),
        "time": t,
        "sites": N,
        "boundary_cells": int(np.sum(b5)),
        "max_final_corners": int(np.max(c5)),
        "max_boundary_final_corners": int(np.max(c5[b5])),
        "cells_with_more_than_48_final_corners": int(np.sum(c5 > 48)),
        "cells_with_more_than_64_final_corners": int(np.sum(c5 > 64)),
        "independent_area_sum": float(np.sum(a5, dtype=np.longdouble)),
        "domain_area": float(DOMAIN_AREA),
        "max_guard_area_relative_difference": float(np.max(guard_rel)),
        "max_native_area_relative_difference": float(np.max(native_rel)),
        "checks": checks,
        "scope": "Accepted native state only; final clipped polygons do not bound intermediate GPU clipping complexity",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))
    if not all(checks.values()):
        raise SystemExit("Bounded checkpoint geometry check failed")


if __name__ == "__main__":
    main()

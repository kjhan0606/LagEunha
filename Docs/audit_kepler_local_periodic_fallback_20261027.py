"""Count periodic images needed by clipped cells beyond a narrow global halo.

This is a read-only exact-orbit geometry control. It does not modify the
native CUDA code or represent valid exterior gas evolution after wrapping.
"""

import argparse
import gc
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi, cKDTree

from audit_kepler_periodic_halo_geometry_20261027 import load_ordered


BOX = 6.0
EPS = 0.01175
CELL = 0.0625


def image_points(xy, margin):
    n = len(xy)
    points = [xy]
    owners = [np.arange(n, dtype=np.int32)]
    sx_values = [np.zeros(n, dtype=np.int8)]
    sy_values = [np.zeros(n, dtype=np.int8)]
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            if sx == sy == 0:
                continue
            mx = (np.ones(n, dtype=bool) if sx == 0 else
                  (xy[:, 0] >= BOX - margin if sx == -1 else xy[:, 0] <= margin))
            my = (np.ones(n, dtype=bool) if sy == 0 else
                  (xy[:, 1] >= BOX - margin if sy == -1 else xy[:, 1] <= margin))
            idx = np.flatnonzero(mx & my)
            points.append(xy[idx] + (sx * BOX, sy * BOX))
            owners.append(idx.astype(np.int32))
            sx_values.append(np.full(len(idx), sx, dtype=np.int8))
            sy_values.append(np.full(len(idx), sy, dtype=np.int8))
    return (np.vstack(points), np.concatenate(owners),
            np.concatenate(sx_values), np.concatenate(sy_values))


def radii(vor, xy, owners):
    result = np.empty(len(owners), dtype=float)
    for k, i in enumerate(owners):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded owned periodic polygon: {i}")
        result[k] = np.max(np.linalg.norm(vor.vertices[region] - xy[i], axis=1))
    if not np.all(np.isfinite(result) & (result > 0)):
        raise ValueError("Invalid Voronoi radius")
    return result


def areas(vor, owners):
    result = np.empty(len(owners), dtype=float)
    for k, i in enumerate(owners):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded audited polygon: {i}")
        vertices = vor.vertices[region]
        result[k] = .5 * abs(np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1)) -
                              np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1)))
    if not np.all(np.isfinite(result) & (result > 0)):
        raise ValueError("Invalid Voronoi area")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("initial", type=Path)
    ap.add_argument("--time", type=float, required=True)
    ap.add_argument("--base-halo-cells", type=int, default=8)
    ap.add_argument("--image-margin", type=float, default=2.5)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if not (args.time >= 0 and 1 <= args.base_halo_cells <= 64 and
            0 < args.image_margin < BOX / 2):
        raise ValueError("Invalid exact-orbit time, halo, or image margin")
    p0, meta = load_ordered(args.initial)
    if float(meta["t"]) != 0 or float(meta["Lx"]) != BOX or float(meta["Ly"]) != BOX:
        raise ValueError("Expected exact t=0 6x6 native Kepler snapshot")
    n = len(p0)
    dx, dy = p0["x"].astype(float) - 3, p0["y"].astype(float) - 3
    r2 = dx * dx + dy * dy
    phase = args.time * (r2 + EPS * EPS)**(-.75)
    xy = np.column_stack((np.mod(3 + dx * np.cos(phase) - dy * np.sin(phase), BOX),
                          np.mod(3 + dx * np.sin(phase) + dy * np.cos(phase), BOX)))
    edge = np.min(np.column_stack((xy[:, 0], xy[:, 1],
                                   BOX - xy[:, 0], BOX - xy[:, 1])), axis=1)
    selected = np.flatnonzero(edge < args.image_margin / 2)
    base_width = args.base_halo_cells * CELL
    narrow_xy, _, _, _ = image_points(xy, base_width)
    narrow_vor = Voronoi(narrow_xy, qhull_options="Qbb Qc Qz")
    partial_radius = radii(narrow_vor, xy, selected)
    clipped = 2 * partial_radius >= edge[selected] + base_width
    risk = selected[clipped]
    risk_partial_radius = partial_radius[clipped]
    del narrow_vor, narrow_xy, partial_radius
    gc.collect()

    full_xy, owner, sx, sy = image_points(xy, args.image_margin)
    full_vor = Voronoi(full_xy, qhull_options="Qbb Qc Qz")
    full_radius = radii(full_vor, xy, risk)
    max_full = float(np.max(full_radius)) if len(risk) else 0.0
    if 2 * max_full >= args.image_margin:
        raise ValueError("Full periodic image margin does not certify risk polygons")
    boxlow, boxhigh = -base_width, BOX + base_width
    in_base_halo = ((sx == 0) & (sy == 0)) | (
        (full_xy[:, 0] >= boxlow) & (full_xy[:, 0] <= boxhigh) &
        (full_xy[:, 1] >= boxlow) & (full_xy[:, 1] <= boxhigh))
    wide_width = 28 * CELL
    in_wide_halo = ((sx == 0) & (sy == 0)) | (
        (full_xy[:, 0] >= -wide_width) & (full_xy[:, 0] <= BOX + wide_width) &
        (full_xy[:, 1] >= -wide_width) & (full_xy[:, 1] <= BOX + wide_width))
    halo32_width = 32 * CELL
    in_32_halo = ((sx == 0) & (sy == 0)) | (
        (full_xy[:, 0] >= -halo32_width) & (full_xy[:, 0] <= BOX + halo32_width) &
        (full_xy[:, 1] >= -halo32_width) & (full_xy[:, 1] <= BOX + halo32_width))

    tree = cKDTree(full_xy)
    local_extra = set()
    per_cell_extra = []
    for i, radius in zip(risk, risk_partial_radius):
        candidates = tree.query_ball_point(xy[i], 2 * radius * (1 + 1e-12))
        extra = [j for j in candidates if not in_base_halo[j]]
        local_extra.update(extra)
        per_cell_extra.append(len(extra))

    risk_mask = np.zeros(len(full_xy), dtype=bool)
    risk_mask[risk] = True
    ridges = full_vor.ridge_points
    left = ridges[risk_mask[ridges[:, 0]], 1]
    right = ridges[risk_mask[ridges[:, 1]], 0]
    actual_new_support = set(np.concatenate((left, right)).tolist()) - set(
        np.flatnonzero(in_base_halo).tolist())
    if not actual_new_support.issubset(local_extra):
        raise ValueError("Local candidate bound missed a full Voronoi face support")
    full_area = areas(full_vor, risk)
    reduced_indices = np.concatenate((np.flatnonzero(in_base_halo),
                                      np.array(sorted(local_extra), dtype=np.int64)))
    if not np.array_equal(reduced_indices[:n], np.arange(n)):
        raise ValueError("Reduced image set does not preserve original owner indices")
    reduced_vor = Voronoi(full_xy[reduced_indices], qhull_options="Qbb Qc Qz")
    reduced_area = areas(reduced_vor, risk)
    area_error = float(np.max(abs(reduced_area / full_area - 1))) if len(risk) else 0.0
    if area_error >= 1e-8:
        raise ValueError("Selective periodic images did not reconstruct full cell areas")
    max_partial = float(np.max(risk_partial_radius)) if len(risk) else 0.0
    result = {
        "initial": str(args.initial), "time": args.time, "real_sites": n,
        "base_halo_cells": args.base_halo_cells,
        "full_image_margin": args.image_margin,
        "selected_edge_cells": int(len(selected)),
        "clipped_cells_from_narrow_voronoi": int(len(risk)),
        "minimum_initial_radius_of_clipped_cells": float(np.min(np.sqrt(r2[risk]))) if len(risk) else None,
        "max_partial_radius": max_partial,
        "max_full_radius_of_clipped_cells": max_full,
        "full_margin_over_twice_max_full_radius": args.image_margin / (2 * max_full) if max_full else None,
        "base_periodic_image_count": int(np.sum(in_base_halo)) - n,
        "global_28_periodic_image_count": int(np.sum(in_wide_halo)) - n,
        "global_32_periodic_image_count": int(np.sum(in_32_halo)) - n,
        "local_extra_image_union_count": len(local_extra),
        "local_extra_candidates_per_clipped_cell_min_median_max":
            [int(np.min(per_cell_extra)), float(np.median(per_cell_extra)),
             int(np.max(per_cell_extra))] if per_cell_extra else [0, 0, 0],
        "actual_new_voronoi_face_support_image_count": len(actual_new_support),
        "actual_new_supports_covered_by_local_candidates": True,
        "reduced_image_count": int(len(reduced_indices) - n),
        "reduced_vs_full_risk_area_relative_max": area_error,
        "reduced_risk_area_passes_1e_minus_8": True,
        "partial_radius_search_covered_by_image_margin": bool(2 * max_partial < args.image_margin),
        "scope": "Exact-orbit wrapped periodic geometry only. Local candidate union is an offline upper-bound strategy using incomplete-cell radius; it is not a native GPU implementation or valid exterior gas evolution."
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

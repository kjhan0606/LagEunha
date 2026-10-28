"""Independently test whether a periodic Voronoi cell needs a wider GPU halo.

Read only a complete native state. Add copies of generators near the 6x6 box
seams and calculate unweighted SciPy Voronoi polygons for original edge cells.
This diagnoses the geometric halo bound; it does not replay the failed RK4
stage or reproduce the GPU candidate algorithm.
"""

import argparse
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import Voronoi

from audit_kepler_polar_checkpoint_20261026 import read_snapshot


def load_ordered(path):
    with path.open("rb") as stream:
        magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
    if magic != -1 or not (0 < n <= 384 * 384) or stride != 576 or stress != 1:
        raise ValueError(f"Unexpected native Kepler snapshot header: {path}")
    if path.stat().st_size != 16 + n * stride + 52:
        raise ValueError(f"Incomplete or overlong native snapshot: {path}")
    p, meta = read_snapshot(path)
    ids = (p["u4if"].astype(np.uint64) &
           np.uint64((1 << 56) - 1)).astype(np.int64)
    if len(p) != n or not np.array_equal(np.sort(ids), np.arange(n)):
        raise ValueError(f"Snapshot cell IDs incomplete: {path}")
    return p[np.argsort(ids)], meta


def shifted_images(xy, width, margin):
    images = [xy]
    count = [len(xy)]
    x, y = xy[:, 0], xy[:, 1]
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            if sx == 0 and sy == 0:
                continue
            mx = np.ones(len(xy), dtype=bool) if sx == 0 else (
                x >= width - margin if sx == -1 else x <= margin)
            my = np.ones(len(xy), dtype=bool) if sy == 0 else (
                y >= width - margin if sy == -1 else y <= margin)
            selected = xy[mx & my]
            images.append(selected + np.array((sx * width, sy * width)))
            count.append(len(selected))
    return np.vstack(images), count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("endpoint", type=Path)
    parser.add_argument("--image-margin", type=float, default=1.)
    parser.add_argument("--halo-cells", type=float, default=8.)
    parser.add_argument("--ballistic-time", type=float)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    p0, m0 = load_ordered(args.initial)
    p, meta = load_ordered(args.endpoint)
    if len(p0) != len(p):
        raise ValueError("Initial and endpoint cell counts differ")
    width = float(m0["Lx"])
    if width != 6. or float(m0["Ly"]) != width or not 0 < args.image_margin < width / 2:
        raise ValueError("Expected 6x6 periodic Kepler state and valid image margin")
    if args.ballistic_time is None:
        xy = np.column_stack((p["x"].astype(float), p["y"].astype(float)))
        result_time = float(meta["t"])
        kind = "native_saved_generator_positions"
        saved_volume_key = "native_saved_volume"
    else:
        if args.ballistic_time < 0:
            raise ValueError("Ballistic control time must be nonnegative")
        x0 = p0["x"].astype(float) - 3.
        y0 = p0["y"].astype(float) - 3.
        r02 = x0 * x0 + y0 * y0
        phase = args.ballistic_time * (r02 + .01175**2)**(-.75)
        xy = np.column_stack((np.mod(3. + x0 * np.cos(phase) - y0 * np.sin(phase), width),
                              np.mod(3. + x0 * np.sin(phase) + y0 * np.cos(phase), width)))
        result_time = args.ballistic_time
        kind = "wrapped_exact_initial_orbits_geometry_only"
        saved_volume_key = "initial_saved_volume"
    if not np.all((xy >= 0.) & (xy < width)):
        raise ValueError("Owned generators are not in the periodic box")
    full, image_count = shifted_images(xy, width, args.image_margin)
    vor = Voronoi(full, qhull_options="Qbb Qc Qz")
    edge = np.min(np.column_stack((xy[:, 0], xy[:, 1],
                                   width - xy[:, 0], width - xy[:, 1])), axis=1)
    selected = np.flatnonzero(edge < args.image_margin / 2.)
    radius = np.empty(len(selected), dtype=float)
    area = np.empty(len(selected), dtype=float)
    for k, i in enumerate(selected):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded edge polygon: {i}")
        vertices = vor.vertices[region]
        radius[k] = np.max(np.linalg.norm(vertices - xy[i], axis=1))
        area[k] = .5 * abs(np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1))
                              - np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1)))
    if not np.all(np.isfinite(radius) & (radius > 0) & np.isfinite(area) & (area > 0)):
        raise ValueError("Invalid independently reconstructed edge geometry")
    halo_width = args.halo_cells * .0625
    padded_edge = edge[selected] + halo_width
    ratio = 2. * radius / padded_edge
    dangerous = ratio >= 1.
    required_halo_width = 2. * radius - edge[selected]
    candidate_halo_cells = (8, 12, 16, 20, 24, 28, 32, 36, 40)
    violations_by_halo_cells = {
        str(cells): int(np.sum(required_halo_width >= cells * .0625))
        for cells in candidate_halo_cells
    }
    r0 = np.hypot(p0["x"].astype(float) - 3., p0["y"].astype(float) - 3.)
    top = []
    for k in np.argsort(ratio)[-20:][::-1]:
        i = selected[k]
        top.append({"id": int(i), "x": float(xy[i, 0]), "y": float(xy[i, 1]),
                    "initial_radius": float(r0[i]),
                    "physical_edge_distance": float(edge[i]),
                    "independent_circumradius": float(radius[k]),
                    "independent_area": float(area[k]),
                    saved_volume_key: float(p["volume"][i]),
                    "twice_radius_over_padded_edge": float(ratio[k])})
    result = {
        "initial": str(args.initial), "endpoint": str(args.endpoint),
        "time": result_time, "kind": kind, "box_width": width,
        "image_margin": args.image_margin, "image_count_by_shift": image_count,
        "points_in_voronoi": int(len(full)),
        "halo_cells": args.halo_cells, "halo_width": halo_width,
        "edge_cells_audited": int(len(selected)),
        "independent_halo_bound_violations": int(np.sum(dangerous)),
        "maximum_twice_radius_over_padded_edge": float(np.max(ratio)),
        "maximum_required_halo_width_for_audited_cells": float(np.max(required_halo_width)),
        "minimum_integer_halo_cells_for_audited_cells": int(np.floor(np.max(required_halo_width) / .0625) + 1),
        "violations_by_halo_cells_for_audited_cells": violations_by_halo_cells,
        "periodic_image_margin_over_twice_max_radius": float(args.image_margin / (2. * np.max(radius))),
        "radius_p99": float(np.quantile(radius, .99)),
        "radius_max": float(np.max(radius)),
        "top_risk_cells": top,
        "scope": ("Independent unweighted periodic Voronoi geometry at saved accepted state only; failed RK4 k4 stage and GPU stencil not replayed"
                  if args.ballistic_time is None else
                  "Wrapped exact initial circular generator orbits only, including nonphysical continuation across the nonsmooth periodic seam; no native gas evolution or valid exterior solution"),
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in (
        "time", "points_in_voronoi", "edge_cells_audited",
        "independent_halo_bound_violations", "maximum_twice_radius_over_padded_edge",
        "maximum_required_halo_width_for_audited_cells",
        "minimum_integer_halo_cells_for_audited_cells",
        "violations_by_halo_cells_for_audited_cells",
        "periodic_image_margin_over_twice_max_radius",
        "radius_p99", "radius_max")}, indent=2))


if __name__ == "__main__":
    main()

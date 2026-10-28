"""Count exact periodic Voronoi corners in an accepted sparse checkpoint."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import IC, digest, snapshot_time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("checkpoint", type=Path)
    ap.add_argument("--expected-count", type=int, choices=tuple(IC), required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    p, meta = ordered(args.checkpoint, args.expected_count)
    t = snapshot_time(args.checkpoint, args.expected_count)
    if abs(meta["t"] - t) > 1e-10:
        raise ValueError("Header/footer time disagreement")
    real = np.column_stack((p["x"].astype(float), p["y"].astype(float)))
    margin = 2.
    pieces = [real]
    image_counts = {}
    for sx in (-6., 0., 6.):
        for sy in (-6., 0., 6.):
            if sx == sy == 0.:
                continue
            shifted = real + (sx, sy)
            take = np.all((shifted > -margin) & (shifted < 6. + margin), axis=1)
            pieces.append(shifted[take])
            image_counts[f"{sx:g},{sy:g}"] = int(np.sum(take))
    points = np.concatenate(pieces)
    vor = Voronoi(points, qhull_options="Qbb Qc Qz")
    corners = np.empty(len(real), dtype=int)
    areas = np.empty(len(real), dtype=float)
    for i in range(len(real)):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded real cell {i} with a two-unit halo")
        corners[i] = len(region)
        vertex = vor.vertices[region]
        areas[i] = .5 * abs(np.dot(vertex[:, 0], np.roll(vertex[:, 1], -1)) -
                              np.dot(vertex[:, 1], np.roll(vertex[:, 0], -1)))
    if not np.all(np.isfinite(areas) & (areas > 0)):
        raise ValueError("Nonpositive or nonfinite periodic Voronoi area")
    volume = p["volume"].astype(float)
    relative_area = np.abs(areas / volume - 1.)
    largest = np.argsort(corners)[-12:][::-1]
    report = {
        "checkpoint": str(args.checkpoint.resolve()),
        "checkpoint_sha256": digest(args.checkpoint),
        "time": t, "real_count": len(real),
        "periodic_image_count": int(len(points) - len(real)),
        "periodic_image_counts_by_shift": image_counts,
        "corners_max": int(np.max(corners)),
        "corners_p99": float(np.quantile(corners, .99)),
        "corners_ge_32": int(np.sum(corners >= 32)),
        "corners_ge_48": int(np.sum(corners >= 48)),
        "corners_ge_64": int(np.sum(corners >= 64)),
        "independent_area_sum": float(np.sum(areas, dtype=np.longdouble)),
        "native_area_sum": float(np.sum(volume, dtype=np.longdouble)),
        "area_relative_max": float(np.max(relative_area)),
        "largest_cells": [{"id": int(i), "corners": int(corners[i]),
                           "x": float(real[i, 0]), "y": float(real[i, 1]),
                           "r": float(np.hypot(real[i, 0]-3., real[i, 1]-3.)),
                           "area_relative": float(relative_area[i])}
                          for i in largest],
        "scope": "Independent Qhull periodic Voronoi at accepted checkpoint; does not reconstruct failed RK4 stage or intermediate clipping polygons",
    }
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps({k: report[k] for k in ("time", "real_count", "corners_max",
                                            "corners_ge_48", "area_relative_max")}))


if __name__ == "__main__":
    main()

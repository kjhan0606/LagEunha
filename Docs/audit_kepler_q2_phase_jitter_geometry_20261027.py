"""Test how the q2 Kepler geometry responds to small angular IC jitter.

This follows exact circular point orbits and recomputes unweighted Voronoi
areas. It is a geometry diagnostic, not a hydrodynamic or heating test.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi

from audit_kepler_polar_checkpoint_20261026 import CX, CY, EPS, load_ordered


TIMES = (0., 6.5, 600.)
BINS = ((.5, .75), (.75, 1.), (1., 2.), (.5, 2.))


def disk_areas(xy, disk):
    vor = Voronoi(xy, qhull_options="Qbb Qc Qz")
    areas = np.empty(len(disk), dtype=float)
    safety = np.empty(len(disk), dtype=float)
    for k, i in enumerate(disk):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded disk Voronoi cell: {i}")
        vertex = vor.vertices[region]
        areas[k] = .5 * abs(np.dot(vertex[:, 0], np.roll(vertex[:, 1], -1))
                              - np.dot(vertex[:, 1], np.roll(vertex[:, 0], -1)))
        circumradius = np.max(np.linalg.norm(vertex - xy[i], axis=1))
        edge = min(xy[i, 0], xy[i, 1], 6. - xy[i, 0], 6. - xy[i, 1])
        safety[k] = edge / (2. * circumradius)
    if not np.all(np.isfinite(areas) & (areas > 0)):
        raise ValueError("Invalid disk cell area")
    if not np.all(safety > 1.):
        raise ValueError("Disk polygons may depend on omitted periodic images")
    return areas, float(np.min(safety))


def one_case(x0, y0, r0, disk, angular_spacing, polar, inverse,
             mode, eta, seed):
    theta0 = np.arctan2(y0, x0)
    if eta:
        rng = np.random.default_rng(seed)
        if mode == "site":
            offset = rng.uniform(-eta, eta, size=len(r0)) * angular_spacing
        else:
            offset = np.zeros(len(r0), dtype=float)
            ring_offset = rng.uniform(-eta, eta, size=int(np.max(inverse))+1)
            offset[polar] = ring_offset[inverse] * angular_spacing[polar]
        theta0 += offset
    initial_area = None
    rows = []
    for t in TIMES:
        phase = theta0 + t * (r0 * r0 + EPS * EPS)**(-.75)
        xy = np.column_stack((CX + r0 * np.cos(phase),
                              CY + r0 * np.sin(phase)))
        area, safety = disk_areas(xy, disk)
        if initial_area is None:
            initial_area = area
        drift = initial_area / area - 1.
        bins = []
        for lo, hi in BINS:
            local = (r0[disk] >= lo) & (r0[disk] < hi)
            weight = initial_area[local]
            bins.append({"rlo": lo, "rhi": hi,
                         "cells": int(np.sum(local)),
                         "density_change_abs_area_weighted": float(
                             np.sum(np.abs(drift[local]) * weight, dtype=np.longdouble)
                             / np.sum(weight, dtype=np.longdouble)),
                         "density_change_p95": float(np.quantile(np.abs(drift[local]), .95)),
                         "density_change_max": float(np.max(np.abs(drift[local])))})
        rows.append({"time": t, "minimum_periodic_image_safety_ratio": safety,
                     "bins": bins})
    return {"jitter_mode": mode,
            "jitter_fraction_of_ring_angular_spacing": eta,
            "seed": seed if eta else None, "times": rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("--jitter-mode", choices=("site", "ring"),
                        default="site")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    p, meta = load_ordered(args.initial)
    if float(meta["t"]) != 0 or not np.all(p["w2"] == 0):
        raise ValueError("Expected complete native q2 W=0 t=0 snapshot")
    x0 = p["x"].astype(float) - CX
    y0 = p["y"].astype(float) - CY
    r0 = np.hypot(x0, y0)
    disk = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    # Only the disk plus immediately adjacent polar guard rings are perturbed.
    polar = (r0 >= .48) & (r0 < 2.02)
    radii, inverse, counts = np.unique(np.round(r0[polar], 8),
                                       return_inverse=True, return_counts=True)
    if len(radii) != 70 or counts.min() < 200:
        raise ValueError("Unexpected q2 disk/guard ring structure")
    angular_spacing = np.zeros(len(r0), dtype=float)
    angular_spacing[polar] = 2. * np.pi / counts[inverse]
    cases = [one_case(x0, y0, r0, disk, angular_spacing, polar, inverse,
                      args.jitter_mode, eta, 20261027)
             for eta in (0., .05, .10)]
    result = {
        "initial": str(args.initial), "disk_cells": int(len(disk)),
        "jittered_radius_range": [0.48, 2.02],
        "jittered_rings": int(len(radii)),
        "jitter_rule": f"{args.jitter_mode} angular offset uniform in ±eta times each ring's 2π/N spacing; radius and exact circular speed remain unchanged; initial disk mass would be 1.01 times each perturbed Voronoi area",
        "cases": cases,
        "scope": "Exact-orbit Voronoi geometry only; no native hydrodynamics, entropy production, or solver-stability prediction",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"disk_cells": result["disk_cells"],
                      "cases": [{"eta": case["jitter_fraction_of_ring_angular_spacing"],
                                 "disk_density_error_by_time": [row["bins"][-1]["density_change_abs_area_weighted"]
                                                                 for row in case["times"]]}
                                for case in cases]}, indent=2))


if __name__ == "__main__":
    main()

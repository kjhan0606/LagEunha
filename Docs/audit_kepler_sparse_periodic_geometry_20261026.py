"""Full periodic Voronoi-area audit of the 54k-site Kepler glass.

Nine tiled copies make every central-box cell finite without relying on the
native GPU tessellator. This is an initial-geometry check, not gas evolution.
"""

import argparse
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import Voronoi


def glass(path):
    with path.open("rb") as stream:
        n, lx, ly = struct.unpack("<idd", stream.read(20))
        x = np.fromfile(stream, dtype="<f8", count=n)
        y = np.fromfile(stream, dtype="<f8", count=n)
    if (not 0 < n < 384*384 or (lx, ly) != (6., 6.) or len(x) != n or len(y) != n or
            path.stat().st_size != 20+16*n):
        raise ValueError("Unexpected glass format")
    xy = np.column_stack((x, y))
    if not np.all((xy >= 0) & (xy < 6)):
        raise ValueError("Generator outside periodic box")
    return xy


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("glass", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    xy = glass(args.glass)
    tiles = [xy]
    for dx in (-6., 0., 6.):
        for dy in (-6., 0., 6.):
            if dx or dy:
                tiles.append(xy + np.array((dx, dy)))
    vor = Voronoi(np.concatenate(tiles), qhull_options="Qbb Qc Qz")
    n = len(xy)
    area = np.empty(n)
    radius = np.empty(n)
    for i in range(n):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded central-cell region i={i}")
        vertices = vor.vertices[region]
        area[i] = .5*abs(np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1))-
                           np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1)))
        radius[i] = np.max(np.linalg.norm(vertices-xy[i], axis=1))
    if not np.all(np.isfinite(area) & (area > 0) & np.isfinite(radius)):
        raise ValueError("Invalid periodic cell geometry")
    r = np.hypot(xy[:, 0]-3., xy[:, 1]-3.)
    profile = np.where(r < .5, (r/.5)**3,
                       np.where(r <= 2., 1., (1.+(r-2.)/.1)**-3))
    density = .01+profile
    mass = density*area
    bins = []
    for lo, hi in ((0., .1), (.1, .48), (.48, .5), (.5, 2.),
                   (2., 2.02), (2.02, 3.), (3., 5.)):
        selected = (r >= lo) & (r < hi)
        if not np.any(selected):
            continue
        bins.append({
            "rlo": lo, "rhi": hi, "cells": int(np.sum(selected)),
            "area_min": float(np.min(area[selected])),
            "area_p50": float(np.median(area[selected])),
            "area_max": float(np.max(area[selected])),
            "circumradius_max": float(np.max(radius[selected])),
            "mass_sum": float(np.sum(mass[selected])),
            "cell_mass_cv": float(np.std(mass[selected])/np.mean(mass[selected])),
        })
    result = {
        "glass": str(args.glass), "n": n, "periodic_tiles": 9,
        "periodic_area_sum": float(np.sum(area, dtype=np.longdouble)),
        "periodic_area_sum_relative_error": float(np.sum(area, dtype=np.longdouble)/36.-1.),
        "area_min": float(np.min(area)), "area_max": float(np.max(area)),
        "maximum_cell_circumradius": float(np.max(radius)),
        "maximum_circumradius_over_nominal_eight_cell_halo":
            float(np.max(radius)/(.0625*8)),
        "estimated_total_mass": float(np.sum(mass, dtype=np.longdouble)),
        "bins": bins,
        "checks": {
            "all_finite_positive": bool(np.all(np.isfinite(area) & (area > 0))),
            "area_partitions_box": bool(abs(np.sum(area, dtype=np.longdouble)/36.-1.) < 1e-8),
            "every_circumradius_inside_nominal_eight_cell_halo":
                bool(np.max(radius) < .0625*8),
        },
        "scope": "Nine-copy periodic initial Voronoi geometry. The nominal eight-grid-cell halo is a simple geometric scale comparison, not proof that the GPU stencil or hydrodynamics is correct.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))
    if not all(result["checks"].values()):
        raise SystemExit("Periodic sparse geometry check failed")


if __name__ == "__main__":
    main()

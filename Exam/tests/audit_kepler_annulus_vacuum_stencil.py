#!/usr/bin/env python3
"""Count a one-layer vacuum-generator stencil for the Hopkins 256^2 annulus.

This is a lattice geometry audit, not a tessellation or hydrodynamics test.
It assumes the retained regular-grid sites have square W=0 Voronoi cells at
initialization and asks how many missing cardinal sites would be needed to
close their four initial faces.
"""

import json
import math
from collections import Counter, defaultdict


def audit(n=256, box=4.0, r_in=0.5, r_out=2.0):
    dx = box / n
    xmin = -0.5 * box
    fluid = set()

    for j in range(n):
        y = xmin + (j + 0.5) * dx
        for i in range(n):
            x = xmin + (i + 0.5) * dx
            r2 = x * x + y * y
            if r_in * r_in <= r2 <= r_out * r_out:
                fluid.add((i, j))

    # For every fluid site, retain the four Cartesian W=0 Voronoi supports.
    # Missing neighbor coordinates become geometry-only vacuum generators.
    incidence = defaultdict(list)
    directions = ((-1, 0), (1, 0), (0, -1), (0, 1))
    for i, j in fluid:
        for di, dj in directions:
            q = (i + di, j + dj)
            if q not in fluid:
                incidence[q].append((di, dj))

    degree_hist = Counter(len(normals) for normals in incidence.values())
    pair_dot_hist = Counter()
    for normals in incidence.values():
        if len(normals) == 2:
            dot = normals[0][0] * normals[1][0] + normals[0][1] * normals[1][1]
            pair_dot_hist[dot] += 1
        elif len(normals) > 2:
            raise AssertionError("more than two fluid faces share a missing site")

    outside_box = sum(
        i < 0 or i >= n or j < 0 or j >= n for i, j in incidence
    )
    cell_area = dx * dx
    guarded_n = n + 2
    guarded_box = box + 2.0 * dx
    guarded_cell_width = guarded_box / guarded_n
    retained_area = len(fluid) * cell_area
    exact_area = math.pi * (r_out * r_out - r_in * r_in)

    result = {
        "resolution": [n, n],
        "lattice_box": [-0.5 * box, 0.5 * box],
        "annulus_radii": [r_in, r_out],
        "fluid_sites": len(fluid),
        "vacuum_geometry_generators": len(incidence),
        "fluid_to_vacuum_faces": sum(k * v for k, v in degree_hist.items()),
        "generator_fluid_incidence_histogram": {
            str(k): degree_hist[k] for k in sorted(degree_hist)
        },
        "shared_generator_normal_dot_histogram": {
            str(k): pair_dot_hist[k] for k in sorted(pair_dot_hist)
        },
        "generators_outside_lattice_box": outside_box,
        "generator_overhead_fraction": len(incidence) / len(fluid),
        "square_cell_area": cell_area,
        "one_cell_guard_resolution": [guarded_n, guarded_n],
        "one_cell_guard_box_width": guarded_box,
        "one_cell_guard_cell_width": guarded_cell_width,
        "one_cell_guard_mean_cell_area": guarded_cell_width * guarded_cell_width,
        "retained_square_cell_area": retained_area,
        "continuum_annulus_area": exact_area,
        "fractional_area_difference": retained_area / exact_area - 1.0,
    }

    if n == 256 and box == 4.0 and r_in == 0.5 and r_out == 2.0:
        assert result["fluid_sites"] == 48240
        assert result["vacuum_geometry_generators"] == 908
        assert result["fluid_to_vacuum_faces"] == 1280
        assert degree_hist == {1: 536, 2: 372}
        assert pair_dot_hist == {0: 372}
        assert outside_box == 88
        assert guarded_n == 258
        assert guarded_box == 4.03125
        assert guarded_cell_width == dx

    return result


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))

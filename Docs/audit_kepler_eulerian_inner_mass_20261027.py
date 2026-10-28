"""Integrate native fixed cell masses inside the Kepler r=0.5 circle.

Each cell is treated as piecewise constant density m/A on its independently
reconstructed W=0 Voronoi polygon. Native V is checked against A, and the
fractional overlap conserves each complete cell's stored mass exactly. This
diagnostic is Eulerian mass at a snapshot, not a time-integrated flux or
subcell hydro model.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi

from plot_kepler_native_density_long_20261026 import complete_ordered


CENTER = np.array([3., 3.])
RADIUS = .5
BOX = 6.


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def segment_disk_signed_area(a, b, radius):
    """Green-theorem contribution of directed edge a->b clipped by a disk."""
    direction = b - a
    aa = float(np.dot(direction, direction))
    if aa == 0:
        return 0.
    bb = 2. * float(np.dot(a, direction))
    cc = float(np.dot(a, a) - radius * radius)
    roots = [0., 1.]
    discriminant = bb * bb - 4. * aa * cc
    if discriminant > 0:
        square = np.sqrt(discriminant)
        roots.extend(t for t in ((-bb - square) / (2. * aa),
                                 (-bb + square) / (2. * aa)) if 0. < t < 1.)
    roots.sort()
    result = 0.
    for left, right in zip(roots[:-1], roots[1:]):
        if right - left < 1e-14:
            continue
        p = a + left * direction
        q = a + right * direction
        middle = .5 * (p + q)
        if float(np.dot(middle, middle)) <= radius * radius:
            result += .5 * cross(p, q)
        else:
            result += .5 * radius * radius * np.arctan2(cross(p, q),
                                                           float(np.dot(p, q)))
    return result


def polygon_disk_area(vertices, radius=RADIUS):
    value = sum(segment_disk_signed_area(a, b, radius)
                for a, b in zip(vertices, np.roll(vertices, -1, axis=0)))
    return abs(value)


def polygon_area(vertices):
    return .5 * abs(np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1)) -
                      np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1)))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def snapshot_inner_mass(path):
    p, meta = complete_ordered(path)
    if not np.all(p["w2"] == 0):
        raise ValueError("W=0 Voronoi snapshots only")
    for field in ("x", "y", "mass", "volume", "den"):
        if not np.all(np.isfinite(p[field])):
            raise ValueError(f"Nonfinite {field} in {path}")
    if np.any(p["mass"] <= 0) or np.any(p["volume"] <= 0):
        raise ValueError("Nonpositive mass or volume")

    xy = np.column_stack((p["x"].astype(float), p["y"].astype(float)))
    relative = xy - CENTER
    # Use the same radius convention as the native radial auditor. Some q1
    # initial generators lie on r=.5 to roundoff; norm(axis=1) classifies
    # three of them differently from hypot despite identical coordinates.
    radii = np.hypot(relative[:, 0], relative[:, 1])
    boundary_near = np.abs(radii - RADIUS) < 1e-12
    nearest_core_radius = float(np.min(radii))
    if nearest_core_radius >= RADIUS:
        raise ValueError("No generator inside inner circle")
    # Any point |x|<=R is within R+rmin of that core generator. A generator
    # with radius >=2R+rmin cannot own an interior point by the triangle
    # inequality. This bounds the complete candidate list without sampling.
    cutoff = 2. * RADIUS + nearest_core_radius
    selected = np.flatnonzero(radii < cutoff)
    vor = Voronoi(xy, qhull_options="Qbb Qc Qz")

    inner_mass = 0.
    boundary_cells = 0
    contributing_cells = 0
    max_area_error = 0.
    min_periodic_safety = float("inf")
    max_overlap_fraction = 0.
    for i in selected:
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded potential inner-circle cell {i}")
        polygon = vor.vertices[region]
        area = polygon_area(polygon)
        if not np.isfinite(area) or area <= 0:
            raise ValueError(f"Invalid Voronoi area for cell {i}")
        circumradius = float(np.max(np.linalg.norm(polygon - xy[i], axis=1)))
        edge = float(np.min((xy[i, 0], xy[i, 1], BOX - xy[i, 0], BOX - xy[i, 1])))
        safety = edge / (2. * circumradius)
        min_periodic_safety = min(min_periodic_safety, safety)
        if safety <= 1.:
            raise ValueError(f"Cell {i} may depend on omitted periodic images")
        area_error = abs(float(p["volume"][i]) / area - 1.)
        max_area_error = max(max_area_error, area_error)
        if area_error > 1e-4:
            raise ValueError(f"Cell {i} native volume differs from Voronoi area: {area_error}")

        overlap = polygon_disk_area(polygon - CENTER)
        if overlap > area * (1. + 1e-8):
            raise ValueError(f"Circle overlap exceeds cell area for {i}")
        fraction = min(1., overlap / area)
        if fraction > 1e-12:
            contributing_cells += 1
            inner_mass += float(p["mass"][i]) * fraction
            max_overlap_fraction = max(max_overlap_fraction, fraction)
            if fraction < 1. - 1e-12:
                boundary_cells += 1

    density_error = float(np.max(np.abs(p["den"].astype(float) /
                                       (p["mass"].astype(float) / p["volume"].astype(float)) - 1.)))
    if density_error > 1e-6:
        raise ValueError(f"Stored density differs from mass/volume: {density_error}")
    return {
        "snapshot": str(path), "sha256": sha256(path), "time": float(meta["t"]),
        "site_count": len(p), "candidate_sites": len(selected),
        "candidate_cutoff_radius": cutoff, "nearest_core_radius": nearest_core_radius,
        "eulerian_piecewise_constant_mass_inside_r0p5": inner_mass,
        "generator_mass_inside_r0p5": float(np.sum(p["mass"][radii < RADIUS], dtype=np.longdouble)),
        "generator_count_within_1e_minus_12_of_r0p5": int(np.sum(boundary_near)),
        "generator_mass_within_1e_minus_12_of_r0p5":
            float(np.sum(p["mass"][boundary_near], dtype=np.longdouble)),
        "contributing_cells": contributing_cells, "boundary_cells": boundary_cells,
        "max_overlap_fraction": max_overlap_fraction,
        "max_native_to_voronoi_area_relative": max_area_error,
        "minimum_periodic_image_safety_ratio": min_periodic_safety,
        "max_density_mass_volume_relative": density_error,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("endpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    start = snapshot_inner_mass(args.initial)
    end = snapshot_inner_mass(args.endpoint)
    if start["time"] != 0 or end["time"] <= 0 or start["site_count"] != end["site_count"]:
        raise ValueError("Expected same-size initial and evolved native snapshots")
    result = {
        "initial": start, "endpoint": end,
        "eulerian_inner_mass_delta":
            end["eulerian_piecewise_constant_mass_inside_r0p5"] -
            start["eulerian_piecewise_constant_mass_inside_r0p5"],
        "generator_inner_mass_delta":
            end["generator_mass_inside_r0p5"] - start["generator_mass_inside_r0p5"],
        "scope": "Independent W=0 Voronoi circle overlap times fixed native cell mass; piecewise constant Eulerian reconstruction. Does not measure integrated face flux or resolve subcell density gradients.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

"""Make a 5% independent angular site-jitter control of the bounded q2 IC."""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import cKDTree

from make_kepler_bounded_polygon_ic_20261027 import (
    BOUNDARY_RADIUS, BOUNDARY_SIDES, BOX, CENTER, GUARD_RADII,
    bounded_areas,
)


SOURCE_SHA256 = "37fcc5f1600bca9083b90aa2a24f98f973acff517266fe623d07e43fb92558c2"
N = 55351
ETA = 0.05
SEED = 20261027


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path)
    ap.add_argument("output_prefix", type=Path)
    args = ap.parse_args()
    glass_path = Path(str(args.output_prefix) + ".glass.bin")
    area_path = Path(str(args.output_prefix) + ".areas.bin")
    manifest_path = Path(str(args.output_prefix) + ".json")
    if any(p.exists() for p in (glass_path, area_path, manifest_path)):
        raise FileExistsError("Refusing to overwrite bounded jitter IC")
    raw = args.source.read_bytes()
    if digest(raw) != SOURCE_SHA256 or len(raw) != 20 + 16 * N:
        raise ValueError("Source is not the accepted 55,351-site bounded q2 glass")
    if struct.unpack_from("<idd", raw) != (N, BOX, BOX):
        raise ValueError("Unexpected native glass header")
    x = np.frombuffer(raw, dtype="<f8", count=N, offset=20).copy()
    y = np.frombuffer(raw, dtype="<f8", count=N, offset=20 + 8 * N).copy()
    initial_x, initial_y = x.copy(), y.copy()
    rx, ry = x - CENTER[0], y - CENTER[1]
    radius = np.hypot(rx, ry)
    polar = (radius >= .48) & (radius < 2.02)
    disk = (radius >= .5) & (radius < 2.)
    rings, inverse, population = np.unique(
        np.round(radius[polar], 8), return_inverse=True, return_counts=True)
    if len(rings) != 70 or int(population.min()) < 200 or int(np.sum(disk)) != 49038:
        raise ValueError("Unexpected q2 disk ring structure")
    angular_spacing = np.zeros(N)
    angular_spacing[polar] = 2 * np.pi / population[inverse]
    angle_offset = np.random.default_rng(SEED).uniform(-ETA, ETA, N)
    angle_offset *= angular_spacing
    angle = np.arctan2(ry[polar], rx[polar]) + angle_offset[polar]
    x[polar] = CENTER[0] + radius[polar] * np.cos(angle)
    y[polar] = CENTER[1] + radius[polar] * np.sin(angle)
    if not (np.array_equal(x[~polar], initial_x[~polar])
            and np.array_equal(y[~polar], initial_y[~polar])):
        raise ValueError("Outer/core sites changed")
    radius_error = float(np.max(np.abs(np.hypot(x - CENTER[0], y - CENTER[1])
                                       - radius)))
    if radius_error > 2e-15 or np.max(radius) >= BOUNDARY_RADIUS * np.cos(
            np.pi / BOUNDARY_SIDES):
        raise ValueError("Jitter changed radius or placed a site beyond wall")
    xy = np.column_stack((x, y))
    nearest = cKDTree(xy).query(xy, k=2, workers=1)[0][:, 1]
    if float(nearest.min()) <= .005:
        raise ValueError("Jitter created a near-collision")
    areas4, _, clipped4 = bounded_areas(xy, GUARD_RADII[0])
    areas5, _, clipped5 = bounded_areas(xy, GUARD_RADII[1])
    area_rel = float(np.max(np.abs(areas4 / areas5 - 1)))
    polygon_area = (0.5 * BOUNDARY_SIDES * BOUNDARY_RADIUS**2
                    * np.sin(2 * np.pi / BOUNDARY_SIDES))
    area_sum = float(np.sum(areas5, dtype=np.longdouble))
    if area_rel > 1e-8 or abs(area_sum - polygon_area) > 1e-7:
        raise ValueError("Independent bounded Voronoi area check failed")
    glass = (struct.pack("<idd", N, BOX, BOX)
             + x.astype("<f8").tobytes() + y.astype("<f8").tobytes())
    sidecar = struct.pack("<idd", N, BOX, BOX) + areas5.astype("<f8").tobytes()
    report = {
        "source": str(args.source), "source_sha256": digest(raw),
        "glass": str(glass_path), "glass_sha256": digest(glass),
        "areas": str(area_path), "areas_sha256": digest(sidecar),
        "sites": N, "disk_sites": int(np.sum(disk)),
        "jittered_sites": int(np.sum(polar)), "jittered_rings": len(rings),
        "jitter_mode": "independent site-angle offset",
        "jitter_fraction_of_ring_spacing": ETA, "seed": SEED,
        "max_radius_error": radius_error,
        "minimum_pair_distance": float(nearest.min()),
        "boundary_cells_per_guard": [clipped4, clipped5],
        "area_sum": area_sum, "polygon_area": polygon_area,
        "max_guard_area_relative_difference": area_rel,
        "scope": "Same bounded exterior and wall; perturb disk/guard site angles only; recompute W0 one-time fixed masses",
    }
    with glass_path.open("xb") as stream:
        stream.write(glass)
    with area_path.open("xb") as stream:
        stream.write(sidecar)
    with manifest_path.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

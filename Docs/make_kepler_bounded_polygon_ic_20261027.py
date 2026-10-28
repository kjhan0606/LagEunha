"""Prepare a fixed-mass q2 Kepler IC inside a nearly circular fixed domain.

The native disk/core/guard generators are kept byte-for-byte. Only sites
outside r=2.65 are dropped and IDs are compacted in original order. Initial
W=0 masses will be rho_target(x_i) times the independently clipped bounded
Voronoi area. This produces an IC sidecar; it does not enable a native
bounded-domain solver or claim late-time stability.
"""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import Voronoi, cKDTree


SOURCE_SHA256 = "de50b97639a243f9384fc7b4bfc402f9eccb52386f30ecdd3450f7b0acaa9e8f"
BOX = 6.0
CENTER = np.array([3.0, 3.0])
KEEP_RADIUS = 2.65
BOUNDARY_RADIUS = 2.75
BOUNDARY_SIDES = 1024
GUARD_SIDES = 256
GUARD_RADII = (4.0, 5.0)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def polygon_area(poly):
    return float(0.5 * abs(np.dot(poly[:, 0], np.roll(poly[:, 1], -1))
                           - np.dot(poly[:, 1], np.roll(poly[:, 0], -1))))


def clip_halfplane(poly, normal, threshold):
    d = poly @ normal - threshold
    if np.all(d <= 0.0):
        return poly
    out = []
    for k, a in enumerate(poly):
        next_k = (k + 1) % len(poly)
        b = poly[next_k]
        da, db = d[k], d[next_k]
        a_in, b_in = da <= 0.0, db <= 0.0
        if a_in:
            out.append(a)
        if a_in != b_in:
            out.append(a + da / (da - db) * (b - a))
    if len(out) < 3:
        raise ValueError("A physical Voronoi cell vanished at the boundary")
    result = np.asarray(out)
    if np.array_equal(result[0], result[-1]):
        result = result[:-1]
    return result


def bounded_areas(xy, guard_radius):
    theta_guard = np.arange(GUARD_SIDES) * (2 * np.pi / GUARD_SIDES)
    guards = CENTER + guard_radius * np.column_stack(
        (np.cos(theta_guard), np.sin(theta_guard)))
    points = np.vstack((xy, guards))
    vor = Voronoi(points, qhull_options="Qbb Qc Qz")
    theta_boundary = (np.arange(BOUNDARY_SIDES) + 0.5) * (
        2 * np.pi / BOUNDARY_SIDES)
    normals = np.column_stack((np.cos(theta_boundary),
                               np.sin(theta_boundary)))
    apothem = BOUNDARY_RADIUS * np.cos(np.pi / BOUNDARY_SIDES)
    thresholds = normals @ CENTER + apothem
    areas = np.empty(len(xy))
    radii = np.empty(len(xy))
    clipped = 0
    for i in range(len(xy)):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded original site {i} with guard radius {guard_radius}")
        poly = vor.vertices[region]
        if np.max(np.linalg.norm(poly - CENTER, axis=1)) > apothem:
            clipped += 1
            for normal, threshold in zip(normals, thresholds):
                poly = clip_halfplane(poly, normal, threshold)
        areas[i] = polygon_area(poly)
        radii[i] = np.max(np.linalg.norm(poly - xy[i], axis=1))
    if not np.all(np.isfinite(areas) & (areas > 0)):
        raise ValueError("Invalid bounded Voronoi areas")
    return areas, radii, clipped


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output_prefix", type=Path)
    args = parser.parse_args()
    glass_path = Path(str(args.output_prefix) + ".glass.bin")
    area_path = Path(str(args.output_prefix) + ".areas.bin")
    manifest_path = Path(str(args.output_prefix) + ".json")
    if any(path.exists() for path in (glass_path, area_path, manifest_path)):
        raise FileExistsError("Refusing to overwrite a bounded IC artifact")

    source = args.source.read_bytes()
    source_sha = digest(source)
    if source_sha != SOURCE_SHA256:
        raise ValueError("Input is not the accepted q2 sparse-wide glass")
    n, lx, ly = struct.unpack_from("<idd", source)
    if (lx, ly) != (BOX, BOX) or len(source) != 20 + 16 * n:
        raise ValueError("Unexpected native glass header or length")
    x = np.frombuffer(source, dtype="<f8", count=n, offset=20)
    y = np.frombuffer(source, dtype="<f8", count=n, offset=20 + 8 * n)
    r = np.hypot(x - CENTER[0], y - CENTER[1])
    keep = r < KEEP_RADIUS
    xy = np.column_stack((x[keep], y[keep]))
    disk_count = int(np.sum((r[keep] >= .5) & (r[keep] < 2.0)))
    if disk_count != 49038 or len(xy) != 55351:
        raise ValueError("Accepted disk or bounded population changed")
    if np.max(r[keep]) >= BOUNDARY_RADIUS * np.cos(np.pi / BOUNDARY_SIDES):
        raise ValueError("A generator starts outside the polygon inradius")

    estimates = [bounded_areas(xy, guard) for guard in GUARD_RADII]
    areas = estimates[-1][0]
    relative_difference = float(np.max(np.abs(estimates[0][0] / areas - 1.0)))
    area_sum = float(np.sum(areas, dtype=np.longdouble))
    polygon_area_exact = (0.5 * BOUNDARY_SIDES * BOUNDARY_RADIUS**2 *
                          np.sin(2 * np.pi / BOUNDARY_SIDES))
    if abs(area_sum - polygon_area_exact) > 1e-7 or relative_difference > 1e-8:
        raise ValueError("Bounded area partition or independent guard-radius comparison failed")
    nearest = cKDTree(xy).query(xy, k=2, workers=1)[0][:, 1]
    if float(nearest.min()) <= 0:
        raise ValueError("Duplicate bounded IC generators")

    glass = (struct.pack("<idd", len(xy), BOX, BOX) +
             xy[:, 0].astype("<f8").tobytes() +
             xy[:, 1].astype("<f8").tobytes())
    sidecar = (struct.pack("<idd", len(xy), BOX, BOX) +
               areas.astype("<f8").tobytes())
    manifest = {
        "source_glass": str(args.source), "source_sha256": source_sha,
        "glass": str(glass_path), "glass_sha256": digest(glass),
        "area_sidecar": str(area_path), "area_sha256": digest(sidecar),
        "sites": len(xy), "disk_sites": disk_count,
        "core_sites": int(np.sum(r[keep] < .5)),
        "guard_and_exterior_sites": int(np.sum((r[keep] >= 2.) & (r[keep] < KEEP_RADIUS))),
        "max_generator_radius": float(np.max(r[keep])),
        "minimum_pair_separation": float(nearest.min()),
        "box": BOX, "domain_center": CENTER.tolist(),
        "domain_circumradius": BOUNDARY_RADIUS,
        "domain_sides": BOUNDARY_SIDES,
        "domain_area": polygon_area_exact,
        "area_sum": area_sum,
        "area_min": float(np.min(areas)), "area_max": float(np.max(areas)),
        "guard_radii": list(GUARD_RADII),
        "clipped_cells_per_guard_radius": [int(row[2]) for row in estimates],
        "max_relative_area_difference_between_guards": relative_difference,
        "max_cell_circumradius": float(np.max(estimates[-1][1])),
        "format": "native little-endian glass/area sidecar; same 6x6 coordinate frame, contiguous IDs",
        "scope": "Bounded W0 one-time fixed-mass Kepler IC only; native CUDA bounded geometry and boundary work remain to be implemented",
    }
    with glass_path.open("xb") as stream:
        stream.write(glass)
    with area_path.open("xb") as stream:
        stream.write(sidecar)
    with manifest_path.open("x") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

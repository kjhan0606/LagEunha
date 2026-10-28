"""Compute exact initial W=0 periodic Voronoi areas for a Kepler glass.

The native CPU KH initial-volume path has a small local stencil and cannot
assign masses to sparse irregular exterior sites. This is a one-time IC
construction tool, not a replacement for native RK4 tessellation.
"""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import Voronoi

from audit_kepler_local_periodic_fallback_20261027 import image_points


BOX = 6.0
MARGINS = (1.5, 2.5)


def one_margin(xy, margin):
    points, _, _, _ = image_points(xy, margin)
    vor = Voronoi(points, qhull_options="Qbb Qc Qz")
    areas = np.empty(len(xy), dtype=np.float64)
    radii = np.empty(len(xy), dtype=np.float64)
    for i in range(len(xy)):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded periodic cell {i} at margin {margin}")
        vertices = vor.vertices[region]
        areas[i] = 0.5 * abs(np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1)) -
                             np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1)))
        radii[i] = np.max(np.linalg.norm(vertices - xy[i], axis=1))
    if not np.all(np.isfinite(areas) & (areas > 0)):
        raise ValueError("Nonpositive or nonfinite periodic Voronoi area")
    if not np.all(2 * radii < margin):
        raise ValueError("Image margin does not cover all face candidates")
    if abs(float(np.sum(areas, dtype=np.longdouble)) - BOX * BOX) > 1e-8:
        raise ValueError("Periodic Voronoi areas do not tile the 6x6 box")
    return areas, radii, len(points) - len(xy)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("glass", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--glass-sha256", required=True)
    args = parser.parse_args()
    manifest = args.output.with_suffix(".json")
    if args.output.exists() or manifest.exists():
        raise FileExistsError("Refusing to overwrite an area sidecar or manifest")
    raw = args.glass.read_bytes()
    source_sha = hashlib.sha256(raw).hexdigest()
    if source_sha != args.glass_sha256:
        raise ValueError("Glass hash does not match the expected IC")
    n, lx, ly = struct.unpack_from("<idd", raw)
    if lx != BOX or ly != BOX or len(raw) != 20 + 16 * n:
        raise ValueError("Invalid native 6x6 glass format")
    x = np.frombuffer(raw, dtype="<f8", count=n, offset=20)
    y = np.frombuffer(raw, dtype="<f8", count=n, offset=20 + 8 * n)
    xy = np.column_stack((x, y))
    if not np.all(np.isfinite(xy) & (xy >= 0) & (xy < BOX)):
        raise ValueError("Invalid glass positions")
    checks = [one_margin(xy, margin) for margin in MARGINS]
    areas = checks[-1][0]
    relative_difference = float(np.max(abs(checks[0][0] / areas - 1)))
    if relative_difference > 1e-8:
        raise ValueError("Independent periodic image margins disagree")

    with args.output.open("xb") as stream:
        stream.write(struct.pack("<idd", n, BOX, BOX))
        stream.write(areas.astype("<f8", copy=False).tobytes())
    output_sha = hashlib.sha256(args.output.read_bytes()).hexdigest()
    report = {"glass": str(args.glass), "glass_sha256": source_sha,
              "area_file": str(args.output), "area_sha256": output_sha,
              "n": n, "box": BOX, "area_sum": float(np.sum(areas, dtype=np.longdouble)),
              "area_min": float(areas.min()), "area_max": float(areas.max()),
              "image_margins": list(MARGINS),
              "periodic_images_per_margin": [int(row[2]) for row in checks],
              "twice_max_circumradius_per_margin": [float(2 * row[1].max()) for row in checks],
              "max_relative_area_difference_between_margins": relative_difference,
              "format": "little-endian int32 count, float64 Lx, float64 Ly, then count float64 areas indexed by native glass site ID",
              "scope": "Exact initial W=0 periodic Voronoi geometry; native solver still needs fixed-mass and initial-volume gates"}
    with manifest.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

"""Audit actual same-ring and adjacent-ring spacing of a bounded Kepler glass."""

import argparse
import hashlib
import json
import struct
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def quantiles(values):
    return {
        key: float(value)
        for key, value in zip(
            ("min", "p10", "median", "p90", "max"),
            np.quantile(values, (0, 0.1, 0.5, 0.9, 1)),
        )
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("glass", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = args.glass.read_bytes()
    n, lx, ly = struct.unpack_from("<idd", data)
    if len(data) != 20 + 16 * n or (lx, ly) != (6.0, 6.0):
        raise ValueError("Unexpected bounded Kepler glass format")
    x = np.frombuffer(data, dtype="<f8", count=n, offset=20)
    y = np.frombuffer(data, dtype="<f8", count=n, offset=20 + 8 * n)
    xy = np.column_stack((x, y))
    radius = np.hypot(x - 3.0, y - 3.0)
    disk = (radius >= 0.5) & (radius < 2.0)
    disk_xy = xy[disk]
    disk_radius = np.round(radius[disk], 8)
    rings, inverse, population = np.unique(
        disk_radius, return_inverse=True, return_counts=True
    )
    if n != 55351 or int(disk.sum()) != 49038 or len(rings) != 68:
        raise ValueError("Not the accepted bounded q2 disk layout")

    same_ring_chord = 2 * rings * np.sin(np.pi / population)
    radial_gap = np.minimum(
        np.r_[np.inf, np.diff(rings)], np.r_[np.diff(rings), np.inf]
    )
    nearest_cross = np.full(len(disk_xy), np.inf)
    for j in range(len(rings) - 1):
        left = np.flatnonzero(inverse == j)
        right = np.flatnonzero(inverse == j + 1)
        left_to_right = cKDTree(disk_xy[right]).query(
            disk_xy[left], k=1, workers=1
        )[0]
        right_to_left = cKDTree(disk_xy[left]).query(
            disk_xy[right], k=1, workers=1
        )[0]
        nearest_cross[left] = np.minimum(nearest_cross[left], left_to_right)
        nearest_cross[right] = np.minimum(nearest_cross[right], right_to_left)

    radial_to_same = radial_gap / same_ring_chord
    cross_to_same = nearest_cross / same_ring_chord[inverse]
    result = {
        "glass": str(args.glass.resolve()),
        "glass_sha256": hashlib.sha256(data).hexdigest(),
        "sites": n,
        "disk_sites": int(disk.sum()),
        "disk_rings": len(rings),
        "disk_radius_range": [float(rings[0]), float(rings[-1])],
        "radial_gap": quantiles(radial_gap),
        "same_ring_chord": quantiles(same_ring_chord),
        "radial_gap_over_same_chord_per_ring": quantiles(radial_to_same),
        "nearest_adjacent_ring_distance_over_same_chord_per_site": quantiles(
            cross_to_same
        ),
        "disk_sites_with_adjacent_ring_nearer_than_same_ring": int(
            np.count_nonzero(cross_to_same < 1)
        ),
        "scope": "Initial generator distances only; not face force, heating, or late-time stability.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

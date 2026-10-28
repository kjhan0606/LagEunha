"""Create the 5% site-angle jitter control from the accepted full q2 glass.

The native Kepler IC path must use EUNHA_KEPLER_MASS_FROM_CELL_VOLUME=1
because this generator changes initial Voronoi volumes. It preserves the
number/order of sites, radii, Cartesian exterior and original native solver.
"""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import cKDTree


EXPECTED_SOURCE_SHA256 = "3ea6edbf56f882d648abea0bef620dac1564f9146194f332f2566af4412ad4a4"
N = 384 * 384
BOX = 6.0
CENTER = 3.0
ETA = 0.05
SEED = 20261027


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix(".json").exists():
        raise FileExistsError("Refusing to overwrite jitter glass or manifest")

    source_bytes = args.source.read_bytes()
    source_sha = hashlib.sha256(source_bytes).hexdigest()
    if source_sha != EXPECTED_SOURCE_SHA256:
        raise ValueError("Source is not the accepted full q2 glass")
    if len(source_bytes) != 20 + 16 * N:
        raise ValueError("Unexpected q2 glass byte count")
    count, lx, ly = struct.unpack_from("<idd", source_bytes)
    if (count, lx, ly) != (N, BOX, BOX):
        raise ValueError("Unexpected q2 glass header")
    x = np.frombuffer(source_bytes, dtype="<f8", count=N, offset=20).copy()
    y = np.frombuffer(source_bytes, dtype="<f8", count=N, offset=20 + 8 * N).copy()
    x0, y0 = x.copy(), y.copy()
    rx, ry = x - CENTER, y - CENTER
    radius = np.hypot(rx, ry)
    polar = (radius >= 0.48) & (radius < 2.02)
    disk = (radius >= 0.5) & (radius < 2.0)
    rings, inverse, populations = np.unique(
        np.round(radius[polar], 8), return_inverse=True, return_counts=True)
    if len(rings) != 70 or int(populations.min()) < 200 or int(disk.sum()) != 49038:
        raise ValueError("Unexpected q2 disk/guard ring structure")

    angular_spacing = np.zeros(N, dtype=np.float64)
    angular_spacing[polar] = 2 * np.pi / populations[inverse]
    # Match audit_kepler_q2_phase_jitter_geometry_20261027.py exactly:
    # one random draw per original site, then multiply by zero outside rings.
    angular_offset = np.random.default_rng(SEED).uniform(-ETA, ETA, N)
    angular_offset *= angular_spacing
    angle = np.arctan2(ry[polar], rx[polar]) + angular_offset[polar]
    x[polar] = CENTER + radius[polar] * np.cos(angle)
    y[polar] = CENTER + radius[polar] * np.sin(angle)
    if not (np.array_equal(x[~polar], x0[~polar]) and
            np.array_equal(y[~polar], y0[~polar])):
        raise ValueError("Nonpolar sites changed")
    radial_error = np.max(np.abs(np.hypot(x - CENTER, y - CENTER) - radius))
    if radial_error > 2e-15 or not np.all((x > 0) & (x < BOX) & (y > 0) & (y < BOX)):
        raise ValueError("Radius or periodic box contract failed")
    nearest = cKDTree(np.column_stack((x, y))).query(
        np.column_stack((x, y)), k=2, workers=1)[0][:, 1]
    min_pair = float(nearest.min())
    if min_pair < 0.35 * BOX / 384:
        raise ValueError("Jitter created a near-collision")

    with args.output.open("xb") as stream:
        stream.write(struct.pack("<idd", N, BOX, BOX))
        stream.write(x.astype("<f8", copy=False).tobytes())
        stream.write(y.astype("<f8", copy=False).tobytes())
    output_sha = hashlib.sha256(args.output.read_bytes()).hexdigest()
    report = {
        "source": str(args.source), "source_sha256": source_sha,
        "output": str(args.output), "output_sha256": output_sha,
        "n": N, "disk_sites": int(disk.sum()),
        "jittered_sites": int(polar.sum()), "jittered_rings": int(len(rings)),
        "jittered_radius_range": [0.48, 2.02],
        "jitter_mode": "independent site angle", "jitter_fraction": ETA,
        "seed": SEED, "max_radius_error": float(radial_error),
        "min_pair_distance": min_pair,
        "requires_native_env": "EUNHA_KEPLER_MASS_FROM_CELL_VOLUME=1",
        "purpose": "Matched native q2 IC sensitivity control; not a separate physical equilibrium",
    }
    with args.output.with_suffix(".json").open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

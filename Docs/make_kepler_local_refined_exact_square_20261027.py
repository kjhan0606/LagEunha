"""Adjust a refined Kepler ring glass to exactly Nx*Ny sites.

The native KH loader only accepts a glass whose site count equals Nx*Ny.
This tool redistributes a small number of additional angular sites over the
existing inner rings, leaving all original exterior coordinates unchanged.
It does not assign masses or weights.
"""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import cKDTree

from plot_kepler_native_density_long_20261026 import complete_ordered


BOX = 6.0
CENTER = 3.0
R_MATCH = 1.0
EXPECTED_SOURCE_COUNT = 175864


def digest_bytes(data):
    return hashlib.sha256(data).hexdigest()


def balanced_increments(counts, extra):
    counts = np.asarray(counts, dtype=np.int64)
    if extra < 0 or counts.size == 0 or counts.sum() <= 0:
        raise ValueError("Cannot distribute requested site-count correction")
    quotas = extra * counts.astype(np.float64) / float(counts.sum())
    increments = np.floor(quotas).astype(np.int64)
    remaining = int(extra - increments.sum())
    if remaining:
        order = np.argsort(-(quotas - increments), kind="stable")
        increments[order[:remaining]] += 1
    if int(increments.sum()) != extra:
        raise AssertionError("Largest-remainder allocation lost sites")
    return increments


def read_glass(path):
    raw = Path(path).read_bytes()
    n, lx, ly = struct.unpack_from("<idd", raw)
    if (lx, ly) != (BOX, BOX) or len(raw) != 20 + 16 * n:
        raise ValueError("Invalid native 6x6 glass file")
    xy = np.column_stack((np.frombuffer(raw, dtype="<f8", count=n, offset=20),
                          np.frombuffer(raw, dtype="<f8", count=n,
                                        offset=20 + 8 * n)))
    if not np.all(np.isfinite(xy) & (xy >= 0) & (xy < BOX)):
        raise ValueError("Invalid glass coordinates")
    return raw, xy


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source_glass", type=Path)
    ap.add_argument("output_glass", type=Path)
    ap.add_argument("--grid-side", type=int, default=420)
    ap.add_argument("--min-separation", type=float, default=0.004)
    args = ap.parse_args()
    target = args.grid_side**2
    source_manifest_path = args.source_glass.with_suffix(".json")
    output_manifest_path = args.output_glass.with_suffix(".json")
    if args.grid_side <= 0 or args.min_separation <= 0:
        raise ValueError("Grid side and minimum separation must be positive")
    if args.output_glass.exists() or output_manifest_path.exists():
        raise FileExistsError("Refusing to overwrite output glass or manifest")

    source_raw, source_xy = read_glass(args.source_glass)
    source_sha = digest_bytes(source_raw)
    source_manifest = json.loads(source_manifest_path.read_text())
    if (len(source_xy) != EXPECTED_SOURCE_COUNT or
            source_manifest.get("output_count") != EXPECTED_SOURCE_COUNT or
            source_manifest.get("output_glass_sha256") != source_sha or
            source_manifest.get("base_time") != 0.0 or
            source_manifest.get("r_match") != R_MATCH):
        raise ValueError("Unexpected locally refined source candidate")
    if target < len(source_xy):
        raise ValueError("Target Nx*Ny count is below source glass count")

    # Independently confirm the source q2 snapshot is W=0 and that its
    # r>=1 coordinates are still exactly the original exterior sites.
    base_path = Path(source_manifest["base_snapshot"])
    base, base_meta = complete_ordered(base_path)
    if (len(base) != 384 * 384 or float(base_meta["t"]) != 0.0 or
            not np.all(base["w2"] == 0)):
        raise ValueError("Source snapshot is not the expected W=0 q2 t=0 state")
    base_xy = np.column_stack((base["x"].astype(np.float64),
                               base["y"].astype(np.float64)))
    base_r = np.hypot(base_xy[:, 0] - CENTER, base_xy[:, 1] - CENTER)
    expected_exterior = base_xy[base_r >= R_MATCH]

    source_r = np.hypot(source_xy[:, 0] - CENTER,
                        source_xy[:, 1] - CENTER)
    source_exterior = source_xy[source_r >= R_MATCH]
    if not np.array_equal(source_exterior, expected_exterior):
        raise ValueError("Source candidate exterior differs from q2 coordinates")

    center_ids = np.flatnonzero(source_r <= 1e-12)
    if len(center_ids) != 1:
        raise ValueError("Expected exactly one central site")
    eligible = np.flatnonzero((source_r > 1e-12) & (source_r < R_MATCH))
    ordered = eligible[np.argsort(source_r[eligible], kind="stable")]
    sorted_r = source_r[ordered]
    # The hypot result can straddle a decimal rounding boundary for points
    # generated on the same ring. Cluster by a tolerance far below the
    # minimum radial spacing instead of rounding radii to fixed decimals.
    split_at = np.flatnonzero(np.diff(sorted_r) > 1e-8) + 1
    source_ring_ids = np.split(ordered, split_at)
    ring_counts = np.asarray([len(ids) for ids in source_ring_ids], dtype=np.int64)
    ring_radii = np.asarray([float(np.median(source_r[ids]))
                             for ids in source_ring_ids], dtype=np.float64)
    if len(ring_radii) != int(source_manifest["refined_ring_count"]):
        raise ValueError("Unexpected refined ring count")

    n_exterior = int(len(source_exterior))
    target_inner = target - n_exterior
    base_inner = 1 + int(ring_counts.sum())
    extra = target_inner - base_inner
    if extra < 0:
        raise ValueError("Target square count cannot preserve all exterior sites")
    increments = balanced_increments(ring_counts, extra)
    final_counts = ring_counts + increments

    inner_rows = [source_xy[center_ids[0]]]
    for j, (radius, count) in enumerate(zip(ring_radii, final_counts), start=1):
        phase = 0.5 * (j % 2)
        theta = 2.0 * np.pi * (np.arange(int(count), dtype=np.float64) + phase) / int(count)
        inner_rows.extend(zip(CENTER + radius * np.cos(theta),
                              CENTER + radius * np.sin(theta)))
    inner_xy = np.asarray(inner_rows, dtype=np.float64)
    xy = np.vstack((inner_xy, source_exterior))
    if len(xy) != target or len(xy) != args.grid_side * args.grid_side:
        raise AssertionError("Corrected glass does not exactly match Nx*Ny")

    # Check periodic separation and unchanged exterior coordinate order.
    tree = cKDTree(xy, boxsize=BOX)
    nearest = tree.query(xy, k=2, workers=1)[0][:, 1]
    if not np.all(np.isfinite(nearest)) or float(nearest.min()) < args.min_separation:
        raise ValueError(f"Refined sites are too close: dmin={nearest.min():.9g}")
    if not np.array_equal(xy[len(inner_xy):], source_exterior):
        raise AssertionError("Exterior q2 coordinates changed")

    with args.output_glass.open("xb") as stream:
        stream.write(struct.pack("<idd", len(xy), BOX, BOX))
        stream.write(xy[:, 0].astype("<f8", copy=False).tobytes())
        stream.write(xy[:, 1].astype("<f8", copy=False).tobytes())
    output_sha = hashlib.sha256(args.output_glass.read_bytes()).hexdigest()
    final_r = np.hypot(xy[:, 0] - CENTER, xy[:, 1] - CENTER)
    bands = [(0.0, 0.5), (0.5, 0.75), (0.75, 0.85),
             (0.85, 0.93), (0.93, 1.0), (1.0, 2.0),
             (2.0, 2.75), (2.75, 5.0)]
    report = {
        "source_glass": str(args.source_glass),
        "source_glass_sha256": source_sha,
        "output_glass": str(args.output_glass),
        "output_glass_sha256": output_sha,
        "source_count": int(len(source_xy)),
        "output_count": int(len(xy)),
        "target_Nx_Ny": [args.grid_side, args.grid_side],
        "output_count_equals_Nx_times_Ny": bool(len(xy) == target),
        "added_inner_sites": int(extra),
        "ring_count": int(len(ring_radii)),
        "ring_angular_site_increments_min_max": [int(increments.min()),
                                                   int(increments.max())],
        "min_periodic_site_separation": float(nearest.min()),
        "median_periodic_site_separation": float(np.median(nearest)),
        "exterior_site_count": n_exterior,
        "exterior_coordinates_bitwise_preserved": True,
        "radial_counts": [{"lo": lo, "hi": hi,
                           "sites": int(np.sum((final_r >= lo) & (final_r < hi)))}
                          for lo, hi in bands],
        "initial_weights": "all zero in the verified source q2 state; native glass stores positions only",
        "mass": "not assigned here; native startup must use exact periodic areas and then preserve each cell mass",
        "scope": "Exact-Nx*Ny W=0 ring glass candidate; no native-evolution claim"
    }
    with output_manifest_path.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

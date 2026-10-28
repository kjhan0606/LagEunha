"""Compare same-initial-ring and adjacent-initial-ring distances in q2's disk."""

import argparse
import json
import struct
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


N_EXPECTED = 384 * 384
STRIDE_EXPECTED = 576
CX = 3.0
CY = 3.0
R_DISK = (0.5, 2.0)

DTYPE = np.dtype({
    "names": ["u4if", "x", "y"],
    "formats": ["<u8", "<f8", "<f8"],
    "offsets": [8, 24, 32],
    "itemsize": STRIDE_EXPECTED,
})


def load(path: Path):
    with path.open("rb") as stream:
        header = stream.read(16)
        if len(header) != 16:
            raise ValueError(f"Truncated native header: {path}")
        magic, n, stride, stress = struct.unpack("<4i", header)
        if (magic, n, stride, stress) != (-1, N_EXPECTED, STRIDE_EXPECTED, 1):
            raise ValueError(f"Unexpected native layout: {(magic, n, stride, stress)}")
        stream.seek(16 + n * stride)
        footer = stream.read(52)
        if len(footer) != 52:
            raise ValueError(f"Missing native footer: {path}")
        time = struct.unpack_from("<d", footer)[0]
    if path.stat().st_size != 16 + n * stride + 52:
        raise ValueError(f"Incomplete native snapshot: {path}")
    data = np.memmap(path, mode="r", dtype=DTYPE, offset=16, shape=(n,))
    ids = np.asarray(data["u4if"], dtype=np.uint64) & np.uint64((1 << 56) - 1)
    if not np.array_equal(np.sort(ids), np.arange(n, dtype=np.uint64)):
        raise ValueError(f"Cell IDs are not unique and complete: {path}")
    return data, ids.astype(np.int64), time


def quantiles(values: np.ndarray) -> dict:
    return {
        "min": float(np.min(values)),
        "p01": float(np.quantile(values, 0.01)),
        "p05": float(np.quantile(values, 0.05)),
        "median": float(np.median(values)),
        "p95": float(np.quantile(values, 0.95)),
        "max": float(np.max(values)),
    }


def analyze(path: Path, ring_label_by_id: np.ndarray, ring_count: int) -> dict:
    data, ids, time = load(path)
    xy = np.column_stack((data["x"], data["y"]))
    radius = np.hypot(xy[:, 0] - CX, xy[:, 1] - CY)
    labels = ring_label_by_id[ids]
    selected = (radius >= R_DISK[0]) & (radius < R_DISK[1]) & (labels >= 0)
    same = np.full(len(ids), np.nan, dtype=float)
    adjacent = np.full(len(ids), np.nan, dtype=float)

    for shell in range(ring_count):
        shell_rows = np.flatnonzero(labels == shell)
        if shell_rows.size < 2:
            raise ValueError(f"Initial shell {shell} has fewer than two cells")
        same_tree = cKDTree(xy[shell_rows])
        same_distances, _ = same_tree.query(xy[shell_rows], k=2)
        same[shell_rows] = same_distances[:, 1]

        adjacent_rows = np.flatnonzero((labels == shell - 1) | (labels == shell + 1))
        if adjacent_rows.size:
            adjacent_tree = cKDTree(xy[adjacent_rows])
            adjacent[shell_rows] = adjacent_tree.query(xy[shell_rows], k=1)[0]

    if not selected.any() or not np.isfinite(same[selected]).all() or not np.isfinite(adjacent[selected]).all():
        raise ValueError(f"Missing same/adjacent-shell distances in disk: {path}")

    ratio = adjacent[selected] / same[selected]
    result = {
        "snapshot": str(path),
        "time": time,
        "disk_cells": int(selected.sum()),
        "same_shell_nearest_distance": quantiles(same[selected]),
        "adjacent_initial_shell_nearest_distance": quantiles(adjacent[selected]),
        "adjacent_over_same_distance_ratio": quantiles(ratio),
        "fraction_adjacent_shell_closer_than_same_shell": float(np.mean(ratio < 1.0)),
        "radial_bins": [],
    }
    for lo, hi in ((0.5, 1.0), (1.0, 1.5), (1.5, 2.0)):
        mask = selected & (radius >= lo) & (radius < hi)
        local_ratio = adjacent[mask] / same[mask]
        result["radial_bins"].append({
            "r": [lo, hi],
            "cells": int(mask.sum()),
            "same_shell_nearest_median": float(np.median(same[mask])),
            "adjacent_shell_nearest_median": float(np.median(adjacent[mask])),
            "adjacent_over_same_median": float(np.median(local_ratio)),
            "fraction_adjacent_closer": float(np.mean(local_ratio < 1.0)),
        })
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("snapshots", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    initial, initial_ids, initial_time = load(args.initial)
    initial_radius = np.hypot(
        np.asarray(initial["x"], dtype=float) - CX,
        np.asarray(initial["y"], dtype=float) - CY,
    )
    in_disk = (initial_radius >= R_DISK[0]) & (initial_radius < R_DISK[1])
    rounded_radius = np.round(initial_radius, 8)
    shell_radii = np.unique(rounded_radius[in_disk])
    labels_at_rows = np.full(N_EXPECTED, -1, dtype=np.int16)
    for shell, ring_radius in enumerate(shell_radii):
        labels_at_rows[rounded_radius == ring_radius] = shell
    label_by_id = np.empty(N_EXPECTED, dtype=np.int16)
    label_by_id[initial_ids] = labels_at_rows
    ring_counts = [int(np.sum(labels_at_rows == shell)) for shell in range(len(shell_radii))]
    ring_spacing = np.diff(shell_radii)

    snapshots = [analyze(path, label_by_id, len(shell_radii)) for path in args.snapshots]
    if any(item["time"] < initial_time for item in snapshots):
        raise ValueError("All analyzed snapshots must be at or after the initial snapshot")
    if any(a["time"] >= b["time"] for a, b in zip(snapshots, snapshots[1:])):
        raise ValueError("Analyzed snapshots must be supplied in strictly increasing time order")

    result = {
        "scope": (
            "Initial polar-ring labels are assigned from radius rounded to 1e-8 "
            "in the t=.001 IC. Distances are Euclidean distances between real "
            "generators; only current centers in 0.5<=r<2 are summarized. This "
            "does not infer fluid face fluxes or prove a t>6 causal mechanism."
        ),
        "initial_snapshot": str(args.initial),
        "initial_time": initial_time,
        "initial_shell_count_in_disk": int(len(shell_radii)),
        "initial_shell_radius_spacing": quantiles(ring_spacing),
        "initial_shell_population": {
            "min": int(min(ring_counts)),
            "median": float(np.median(ring_counts)),
            "max": int(max(ring_counts)),
        },
        "snapshots": snapshots,
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

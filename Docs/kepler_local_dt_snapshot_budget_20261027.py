"""Read-only local-timestep arithmetic audit for a native 384^2 KH snapshot.

The serialized per-cell `dt` is the final-stage hydro CFL candidate at byte
128 of the 576-byte stress-particle record. This script combines it with the
solver's softened orbital bound. Power-of-two update counts are frozen-state
arithmetic only; they do not model conservative face registers, mesh cost,
positivity, or wall time.
"""

import argparse
import json
import struct
from pathlib import Path

import numpy as np


N_EXPECTED = 384 * 384
STRIDE_EXPECTED = 576


def load(path: Path):
    with path.open("rb") as stream:
        header = stream.read(16)
        if len(header) != 16:
            raise ValueError(f"Truncated native header: {path}")
        magic, n, stride, stress = struct.unpack("<4i", header)
        if (magic, n, stride, stress) != (-1, N_EXPECTED, STRIDE_EXPECTED, 1):
            raise ValueError(f"Unexpected snapshot layout: {(magic, n, stride, stress)}")
        stream.seek(16 + n * stride)
        footer = stream.read(8)
        if len(footer) != 8:
            raise ValueError(f"Missing native time footer: {path}")
        time = struct.unpack("<d", footer)[0]

    dtype = np.dtype({
        "names": ["u4if", "x", "y", "dt"],
        "formats": ["<u8", "<f8", "<f8", "<f4"],
        "offsets": [8, 24, 32, 128],
        "itemsize": stride,
    })
    data = np.memmap(path, mode="r", dtype=dtype, offset=16, shape=(n,))
    ids = data["u4if"].astype(np.uint64) & np.uint64((1 << 56) - 1)
    if not np.array_equal(np.sort(ids), np.arange(n, dtype=np.uint64)):
        raise ValueError(f"Cell IDs are not unique and complete: {path}")
    hydro = np.asarray(data["dt"], dtype=np.float64)
    if not (np.isfinite(hydro).all() and (hydro > 0).all()):
        raise ValueError(f"Invalid local hydro timestep: {path}")
    radius = np.hypot(data["x"] - 3.0, data["y"] - 3.0)
    return data, time, radius, hydro


def weighted_stats(values, mask):
    selected = values[mask]
    if selected.size == 0:
        return {"cells": 0, "min": None, "median": None, "p95": None}
    return {
        "cells": int(selected.size),
        "min": float(np.min(selected)),
        "median": float(np.median(selected)),
        "p95": float(np.quantile(selected, 0.95)),
    }


def analyze(snapshot: Path, eta: float, eps: float, spans: list[float]) -> dict:
    data, time, radius, hydro = load(snapshot)
    orbital = eta * (radius * radius + eps * eps) ** 0.75
    local = np.minimum(hydro, orbital)
    base = float(np.min(local))
    levels = np.floor(np.log2(local / base)).astype(np.int32)
    block_dt = base * np.exp2(levels)
    disk = (radius >= 0.5) & (radius < 2.0)
    regions = []
    edges = (0.0, 0.1, 0.5, 0.75, 1.0, 2.0, 2.35, 2.75, np.inf)
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (radius >= lo) & (radius < hi)
        regions.append({
            "r0": [lo, hi if np.isfinite(hi) else None],
            "cells": int(np.sum(mask)),
            "combined_dt_min": float(np.min(local[mask])) if np.any(mask) else None,
            "combined_dt_median": float(np.median(local[mask])) if np.any(mask) else None,
            "hydro_dt_min": float(np.min(hydro[mask])) if np.any(mask) else None,
            "orbital_dt_min": float(np.min(orbital[mask])) if np.any(mask) else None,
        })

    horizon_results = []
    n = len(local)
    for span in spans:
        global_steps = int(np.ceil(span / base))
        block_updates = np.ceil(span / block_dt).astype(np.int64)
        ideal_updates = np.ceil(span / local).astype(np.int64)
        level_hist = []
        for level in np.unique(levels):
            mask = levels == level
            level_hist.append({
                "level": int(level),
                "cells": int(np.sum(mask)),
                "active_cell_updates": int(np.sum(block_updates[mask], dtype=np.int64)),
            })
        horizon_results.append({
            "frozen_span": span,
            "global_steps": global_steps,
            "all_cell_updates_if_global": int(n * global_steps),
            "ideal_local_cell_updates": int(np.sum(ideal_updates, dtype=np.int64)),
            "power2_block_cell_updates": int(np.sum(block_updates, dtype=np.int64)),
            "power2_block_arithmetic_ratio": float(
                n * global_steps / np.sum(block_updates, dtype=np.int64)),
            "levels": level_hist,
        })

    return {
        "scope": (
            "Frozen-state arithmetic from one native checkpoint; no wall-time "
            "or valid asynchronous hydro integrator is implied."
        ),
        "snapshot": str(snapshot),
        "time": time,
        "cells": n,
        "eta": eta,
        "eps": eps,
        "hydro_dt_min": float(np.min(hydro)),
        "hydro_dt_median": float(np.median(hydro)),
        "orbital_dt_min": float(np.min(orbital)),
        "combined_global_dt": base,
        "cells_within_2x_global_dt": int(np.sum(local < 2 * base)),
        "cells_within_4x_global_dt": int(np.sum(local < 4 * base)),
        "disk_0p5_to_2_combined_dt": weighted_stats(local, disk),
        "radial_regions_by_initial_position": regions,
        "horizons": horizon_results,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--eta", type=float, default=0.06)
    parser.add_argument("--eps", type=float, default=0.01175)
    parser.add_argument("--span", type=float, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.eta <= 0 or args.eps <= 0 or any(span <= 0 for span in args.span):
        raise ValueError("eta, eps and every span must be positive")
    result = analyze(args.snapshot, args.eta, args.eps, args.span)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

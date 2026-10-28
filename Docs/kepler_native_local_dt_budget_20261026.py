"""Measure cellwise hydro/orbital timestep headroom in native KH snapshots.

The serialized `dt` field is a float32 at byte 128 of the 576-byte
treevorostressrk4particletype: see eunha.h VoroQ after the double x/y/mass/v/a
fields. The run's final GPU-step minimum is checked against this field.
"""

import argparse
import json
import re
import struct
from pathlib import Path

import numpy as np


def read_native(path: Path) -> tuple[np.memmap, float]:
    with path.open("rb") as stream:
        magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
        if (magic, n, stride, stress) != (-1, 384*384, 576, 1):
            raise ValueError(f"Unexpected snapshot layout: {(magic,n,stride,stress)}")
        stream.seek(16+n*stride)
        t = struct.unpack("<d", stream.read(8))[0]
    dtype = np.dtype({"names": ["u4if", "x", "y", "dt"],
                      "formats": ["<u8", "<f8", "<f8", "<f4"],
                      "offsets": [8, 24, 32, 128], "itemsize": stride})
    data = np.memmap(path, mode="r", dtype=dtype, offset=16, shape=(n,))
    ids = data["u4if"].astype(np.uint64) & np.uint64((1 << 56)-1)
    if not np.array_equal(np.sort(ids), np.arange(n, dtype=np.uint64)):
        raise ValueError(f"Cell IDs are not unique and complete: {path}")
    return data, t


def final_gpu_dt(log: Path) -> float:
    pattern = re.compile(r"\[GPU step \d+\].*\bdt=([0-9.eE+-]+)")
    matches = pattern.findall(log.read_text(errors="replace"))
    if not matches:
        raise ValueError(f"No GPU step dt found in {log}")
    return float(matches[-1])


def analyze(snapshot: Path, eta: float, eps: float, span: float) -> dict:
    p, t = read_native(snapshot)
    hydro = np.asarray(p["dt"], dtype=np.float64)
    if not (np.isfinite(hydro).all() and (hydro > 0).all()):
        raise ValueError(f"Invalid local hydro timestep: {snapshot}")
    logged = final_gpu_dt(snapshot.parent / "log")
    if abs(float(np.min(hydro))/logged-1) > 5e-4:
        raise ValueError(f"Serialized hydro minimum disagrees with GPU log: {snapshot}")
    r = np.hypot(p["x"]-3.0, p["y"]-3.0)
    orbital = eta * (r*r+eps*eps)**0.75
    local = np.minimum(hydro, orbital)
    base = float(np.min(local))
    n = len(local)
    global_steps = int(np.ceil(span/base))
    levels = np.floor(np.log2(local/base)).astype(np.int32)
    block_dt = base * np.exp2(levels)
    block_updates = np.ceil(span/block_dt)
    ideal_updates = np.ceil(span/local)
    disk = (r >= .5) & (r < 2.)
    level_hist = []
    for level in np.unique(levels):
        mask = levels == level
        level_hist.append({"level": int(level), "cells": int(np.sum(mask)),
                           "active_cell_updates": int(np.sum(block_updates[mask]))})
    return {
        "snapshot": str(snapshot), "time": t, "cells": n,
        "hydro_min_from_serialized_cells": float(np.min(hydro)),
        "hydro_min_from_final_gpu_log": logged,
        "orbital_min": float(np.min(orbital)),
        "actual_global_bound_from_two_terms": base,
        "disk_cells": int(np.sum(disk)),
        "disk_local_hydro_dt_min_median": [float(np.min(hydro[disk])),
                                           float(np.median(hydro[disk]))],
        "disk_local_combined_dt_min_median": [float(np.min(local[disk])),
                                              float(np.median(local[disk]))],
        "cells_within_2x_global_bound": int(np.sum(local < 2*base)),
        "cells_within_4x_global_bound": int(np.sum(local < 4*base)),
        "global_steps_for_span_at_snapshot_dt": global_steps,
        "all_cell_updates_if_global": int(n*global_steps),
        "ideal_local_cell_updates": int(np.sum(ideal_updates)),
        "power2_block_cell_updates": int(np.sum(block_updates)),
        "power2_block_arithmetic_ratio": float(n*global_steps/np.sum(block_updates)),
        "levels": level_hist,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("snapshots", nargs="+", type=Path)
    ap.add_argument("--eta", type=float, default=.06)
    ap.add_argument("--eps", type=float, default=.01175)
    ap.add_argument("--span", type=float, default=.825)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    if args.eta <= 0 or args.eps <= 0 or args.span <= 0:
        raise ValueError("eta, eps, span must be positive")
    result = {"scope": "Native t=0.001 K4 hydro dt and current orbital bound, "
               "frozen over a hypothetical t=0.825 span. This is active-cell "
               "arithmetic only, not a measured or guaranteed runtime speedup.",
              "cases": [analyze(s, args.eta, args.eps, args.span)
                        for s in args.snapshots]}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

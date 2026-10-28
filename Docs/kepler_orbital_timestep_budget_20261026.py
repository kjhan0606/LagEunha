"""Read native polar glass files and quantify the global Kepler dt bottleneck.

This is an orbital-limiter work count, not a hydro CFL or runtime prediction.
It deliberately leaves the solver and existing runs untouched.
"""

import argparse
import json
import struct
from pathlib import Path

import numpy as np


def read_glass(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open("rb") as stream:
        n, lx, ly = struct.unpack("<idd", stream.read(20))
        if (n, lx, ly) != (384 * 384, 6.0, 6.0):
            raise ValueError(f"Unexpected glass header in {path}: {(n, lx, ly)}")
        x = np.fromfile(stream, dtype="<f8", count=n)
        y = np.fromfile(stream, dtype="<f8", count=n)
        if len(x) != n or len(y) != n or stream.read(1):
            raise ValueError(f"Incomplete or overlong glass file: {path}")
    return x, y


def summarize(path: Path, eta: float, eps: float, span: float,
              hydro_cap: float | None) -> dict:
    x, y = read_glass(path)
    r = np.hypot(x - 3.0, y - 3.0)
    local_dt = eta * (r * r + eps * eps) ** 0.75
    global_dt = float(np.min(local_dt))
    disk = (r >= 0.5) & (r < 2.0)
    if not np.any(disk):
        raise ValueError("No disk cells")

    # Power-of-two blocks at the same base dt. This estimates only active-cell
    # work; a conservative local face update, mesh build, and global sync are
    # not included. A floor keeps each chosen block within the orbital bound.
    levels = np.floor(np.log2(local_dt / global_dt)).astype(np.int32)
    block_dt = global_dt * np.exp2(levels)
    ideal_updates = np.ceil(span / local_dt)
    block_updates = np.ceil(span / block_dt)
    global_updates = int(np.ceil(span / global_dt))
    n = len(r)
    bands = [(0.0, 0.1), (0.1, 0.5), (0.5, 0.75),
             (0.75, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, 5.0)]
    radial = []
    for lo, hi in bands:
        mask = (r >= lo) & (r < hi)
        if np.any(mask):
            radial.append({"r_lo": lo, "r_hi": hi, "cells": int(np.sum(mask)),
                           "orbital_dt_min": float(np.min(local_dt[mask])),
                           "orbital_dt_median": float(np.median(local_dt[mask])),
                           "block_level_min": int(np.min(levels[mask])),
                           "block_level_median": float(np.median(levels[mask]))})
    result = {
        "glass": str(path), "cells": n, "eta": eta, "eps": eps,
        "span": span, "minimum_radius": float(np.min(r)),
        "global_orbital_dt": global_dt,
        "disk_cells": int(np.sum(disk)),
        "disk_minimum_orbital_dt": float(np.min(local_dt[disk])),
        "disk_minimum_dt_over_global": float(np.min(local_dt[disk]) / global_dt),
        "cells_within_2x_global_orbital_dt": int(np.sum(local_dt < 2 * global_dt)),
        "cells_within_4x_global_orbital_dt": int(np.sum(local_dt < 4 * global_dt)),
        "global_steps_over_span": global_updates,
        "global_cell_updates_over_span": n * global_updates,
        "ideal_local_cell_updates_over_span": int(np.sum(ideal_updates)),
        "power2_block_cell_updates_over_span": int(np.sum(block_updates)),
        "orbital_only_ideal_work_ratio": float(n * global_updates / np.sum(ideal_updates)),
        "orbital_only_power2_work_ratio": float(n * global_updates / np.sum(block_updates)),
        "radial_bands": radial,
    }
    if hydro_cap is not None:
        if hydro_cap <= 0:
            raise ValueError("Hydro cap must be positive")
        capped_dt = np.minimum(local_dt, hydro_cap)
        capped_levels = np.floor(np.log2(capped_dt / global_dt)).astype(np.int32)
        capped_blocks = global_dt * np.exp2(capped_levels)
        capped_ideal_updates = np.ceil(span / capped_dt)
        capped_block_updates = np.ceil(span / capped_blocks)
        result["global_hydro_cfl_observed_initial"] = hydro_cap
        result["orbital_and_global_hydro_cap_ideal_cell_work_ratio"] = float(
            n * global_updates / np.sum(capped_ideal_updates))
        result["orbital_and_global_hydro_cap_power2_cell_work_ratio"] = float(
            n * global_updates / np.sum(capped_block_updates))
        result["orbital_and_global_hydro_cap_power2_cell_updates"] = int(
            np.sum(capped_block_updates))
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("glasses", nargs="+", type=Path)
    ap.add_argument("--eta", type=float, default=0.06)
    ap.add_argument("--eps", type=float, default=0.01175)
    ap.add_argument("--span", type=float, default=0.825)
    ap.add_argument("--hydro-caps", nargs="*", default=None,
                    help="Initial observed global hydro CFL per glass, or 'none'")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    if args.eta <= 0 or args.eps <= 0 or args.span <= 0:
        raise ValueError("eta, eps, and span must be positive")
    if args.hydro_caps is not None and len(args.hydro_caps) != len(args.glasses):
        raise ValueError("Supply one hydro cap or 'none' per glass")
    caps = ([None] * len(args.glasses) if args.hydro_caps is None else
            [None if value.lower() == "none" else float(value)
             for value in args.hydro_caps])
    result = {"scope": "Initial native polar generators. Hydro caps, when given, "
               "are observed first-step global hydro CFL values applied to all cells. "
               "Local hydro CFL, mesh cost, and conservative face synchronization "
               "are not modeled; work ratios are not runtime predictions.",
              "cases": [summarize(p, args.eta, args.eps, args.span, cap)
                        for p, cap in zip(args.glasses, caps)]}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

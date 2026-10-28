"""Locate saved native hydro/CFL limits by initial region and cell ID."""

import argparse
import json
import struct
from pathlib import Path

import numpy as np


DTYPE = np.dtype({"names": ("tag", "x", "y", "dt"),
                  "formats": ("<u8", "<f8", "<f8", "<f4"),
                  "offsets": (8, 24, 32, 128), "itemsize": 576})


def read(path: Path):
    with path.open("rb") as stream:
        magic, count, stride, stress = struct.unpack("<4i", stream.read(16))
        stream.seek(-52, 2)
        time = struct.unpack("<d", stream.read(8))[0]
    if (magic, stride, stress) != (-1, 576, 1) or path.stat().st_size != 16 + 576*count + 52:
        raise ValueError(f"Incomplete native snapshot: {path}")
    raw = np.memmap(path, mode="r", dtype=DTYPE, offset=16, shape=(count,))
    ids = raw["tag"].astype(np.uint64) & np.uint64((1 << 56)-1)
    order = np.argsort(ids)
    if not np.array_equal(ids[order], np.arange(count, dtype=np.uint64)):
        raise ValueError(f"Incomplete or duplicate IDs: {path}")
    return raw[order], time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("initial", type=Path)
    ap.add_argument("checkpoint", type=Path)
    ap.add_argument("--eta", type=float, default=.06)
    ap.add_argument("--eps", type=float, default=.01175)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    p0, t0 = read(args.initial)
    p, t = read(args.checkpoint)
    if len(p0) != len(p) or abs(t0) > 1e-12 or not (t > t0):
        raise ValueError("Initial and checkpoint times/count mismatch")
    dt = p["dt"].astype(np.float64)
    if not np.all(np.isfinite(dt) & (dt > 0)):
        raise ValueError("Nonpositive or nonfinite stored hydro timestep")
    r0 = np.hypot(p0["x"]-3., p0["y"]-3.)
    radius = np.hypot(p["x"]-3., p["y"]-3.)
    edge0 = np.minimum.reduce((p0["x"], p0["y"], 6.-p0["x"], 6.-p0["y"]))
    edge = np.minimum.reduce((p["x"], p["y"], 6.-p["x"], 6.-p["y"]))
    orbital = args.eta*(radius*radius+args.eps*args.eps)**.75
    combined = np.minimum(dt, orbital)
    def row(i):
        return {"id": int(i), "hydro_dt": float(dt[i]),
                "orbital_dt": float(orbital[i]), "combined_dt": float(combined[i]),
                "initial_r": float(r0[i]), "current_r": float(radius[i]),
                "initial_edge_distance": float(edge0[i]),
                "current_edge_distance": float(edge[i])}
    groups = {
        "core_r0_lt_0p5": r0 < .5,
        "disk_r0_0p5_to_2": (r0 >= .5) & (r0 < 2.),
        "guard_r0_2_to_2p35": (r0 >= 2.) & (r0 < 2.35),
        "buffer_r0_2p35_to_2p75": (r0 >= 2.35) & (r0 < 2.75),
        "outer_r0_ge_2p75": r0 >= 2.75,
        "seam_edge0_lt_0p2": edge0 < .2,
    }
    origin = {}
    for name, mask in groups.items():
        sub = np.flatnonzero(mask)
        origin[name] = {"cells": int(len(sub)), "min_hydro_dt": float(np.min(dt[sub])),
                        "min_orbital_dt": float(np.min(orbital[sub])),
                        "within_2x_global_bound": int(np.sum(combined[sub] < 2*np.min(combined)))}
    result = {"initial": str(args.initial), "checkpoint": str(args.checkpoint),
              "time": t, "cells": len(p), "eta": args.eta, "eps": args.eps,
              "min_hydro_dt": float(np.min(dt)),
              "min_orbital_dt": float(np.min(orbital)),
              "min_combined_dt": float(np.min(combined)),
              "min_hydro_cell": row(int(np.argmin(dt))),
              "min_orbital_cell": row(int(np.argmin(orbital))),
              "min_combined_cell": row(int(np.argmin(combined))),
              "smallest_hydro_cells": [row(int(i)) for i in np.argsort(dt)[:12]],
              "groups": origin,
              "scope": "Saved dt is a per-cell final-RK-stage hydro candidate at this checkpoint, not proof of which cell set a preceding transient accepted global dt. The orbital bound is reconstructed from current generator radius."}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in
                      ("time", "min_hydro_dt", "min_orbital_dt", "min_combined_dt",
                       "min_hydro_cell", "min_combined_cell")}, indent=2))


if __name__ == "__main__":
    main()

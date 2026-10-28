"""Read-only late-checkpoint CFL locality; no mesh or trajectory evolution."""

import argparse
from collections import deque
import hashlib
import json
from pathlib import Path
import re
import struct
import sys

import numpy as np

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_checked(path):
    with path.open("rb") as stream:
        magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
    if ((magic, n, stride, stress) != (-1, 384*384, 576, 1)
            or path.stat().st_size != 16+n*stride+52):
        raise ValueError("Incomplete or incompatible native snapshot")
    p, meta = read_snapshot(path)
    ids = (p["u4if"].astype(np.uint64) & np.uint64((1 << 56)-1)).astype(np.int64)
    order = np.argsort(ids)
    if not np.array_equal(ids[order], np.arange(n)):
        raise ValueError("Incomplete or repeated generator IDs")
    # Compiler-reported VoroQ dt float at byte 128, already used by the
    # earlier native local-dt auditor. Independently cross-check with log.
    dtype = np.dtype({"names": ["dt"], "formats": ["<f4"],
                      "offsets": [128], "itemsize": stride})
    view = np.memmap(path, mode="r", dtype=dtype, offset=16, shape=(n,))
    return p[order], np.asarray(view["dt"], dtype=float)[order], meta


def matched_log(path, step):
    gpu = deque(maxlen=4)
    with path.open(errors="replace") as stream:
        for line in stream:
            match = re.search(r"\[GPU step \d+\].*\bdt=([0-9.eE+-]+)", line)
            if match:
                gpu.append(float(match[1]))
            match = re.search(r"Time=\s*(\S+)\s+step=\s*(\d+)\s+dt=\s*(\S+)", line)
            if match and int(match[2]) == step:
                if len(gpu) != 4:
                    raise ValueError("Missing four preceding CUDA force records")
                return {"step": step, "printed_time": float(match[1]),
                        "accepted_dt": float(match[3]), "last_four_GPU_dt": list(gpu)}
    raise ValueError("Snapshot step absent from log")


def analyze(initial, snapshot, log, step, courant, gamma):
    p0, _, m0 = read_checked(initial)
    p, dt, meta = read_checked(snapshot)
    if m0["t"] != 0 or not meta["t"] > 0 or len(p) != len(p0):
        raise ValueError("Expected matching t=0 and evolved checkpoints")
    if not (np.isfinite(dt).all() and (dt > 0).all()):
        raise ValueError("Nonpositive or nonfinite stored dt")
    row = matched_log(log, step)
    base = float(dt.min())
    # These are rounded native log strings, not new numerical tolerances.
    if abs(row["printed_time"]-meta["t"]) > 5.1e-6:
        raise ValueError("Snapshot time disagrees with matching native step")
    if abs(row["accepted_dt"]/meta["dt"]-1) > 5.1e-6:
        raise ValueError("Snapshot accepted dt disagrees with log")
    if abs(base/row["last_four_GPU_dt"][-1]-1) > 5.1e-4:
        raise ValueError("Stored minimum does not match final CUDA stage")
    for field in ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2"):
        if not all(np.isfinite(q[field]).all() for q in (p0, p)):
            raise ValueError("Nonfinite physical state")
    if not (np.array_equal(p0["mass"], p["mass"]) and np.all(p["w2"] == 0)
            and np.all(p0["w2"] == 0) and np.all(p["den"] > 0)
            and np.all(p["pressure"] > 0)
            and all(np.all(q["ie"] > 0) and np.all(q["volume"] > 0) for q in (p0, p))):
        raise ValueError("Expected fixed-mass positive W0 checkpoint")
    r0 = np.hypot(p0["x"]-3., p0["y"]-3.)
    radius = np.hypot(p["x"]-3., p["y"]-3.)
    levels = np.floor(np.log2(dt/base)).astype(int)
    disk = (r0 >= .5) & (r0 < 2.)
    bins = []
    for lo, hi in ((0., .1), (.1, .5), (.5, .75), (.75, 1.),
                   (1., 2.), (2., 2.35), (2.35, 2.75), (2.75, 5.)):
        mask = (r0 >= lo) & (r0 < hi)
        bins.append({"r0_lo": lo, "r0_hi": hi, "cells": int(mask.sum()),
                     "stored_K4_dt_min": float(dt[mask].min()),
                     "stored_K4_dt_median": float(np.median(dt[mask]))})
    box = np.array([meta["Lx"], meta["Ly"]])
    xy = np.column_stack((p["x"], p["y"]))
    velocity = np.column_stack((p["vx"], p["vy"]))
    limiting = np.flatnonzero(dt < 2*base)
    pairs = []
    # Bound the direct-neighbor work even if another checkpoint has many
    # tied limiters. This audit never reconstructs the full mesh.
    selected = limiting[np.argsort(dt[limiting])[:16]]
    for i in selected:
        displacement = xy-xy[i]
        displacement -= box*np.rint(displacement/box)
        distance2 = np.sum(displacement*displacement, axis=1)
        distance2[i] = np.inf
        j = int(np.argmin(distance2))
        d = float(np.sqrt(distance2[j]))
        if not d > 0:
            raise ValueError("Coincident periodic generators")
        dv = velocity[j]-velocity[i]
        closing = float(-np.dot(dv, displacement[j])/d)
        cs = np.sqrt(gamma*np.array([p["pressure"][i]/p["den"][i],
                                     p["pressure"][j]/p["den"][j]], dtype=float))
        pair_dt = float(2*courant*d/(sum(cs)+max(0., closing)))
        pairs.append({"id": int(i), "initial_radius": float(r0[i]),
                      "current_radius": float(radius[i]), "stored_K4_dt": float(dt[i]),
                      "initial_xy": [float(p0["x"][i]), float(p0["y"][i])],
                      "current_xy": xy[i].tolist(), "current_vxvy": velocity[i].tolist(),
                      "specific_Lz": float((xy[i, 0]-3)*velocity[i, 1]
                                            -(xy[i, 1]-3)*velocity[i, 0]),
                      "nearest_periodic_id": j, "minimum_image_distance": d,
                      "relative_speed": float(np.linalg.norm(dv)),
                      "normal_closing_speed": closing, "sound_speed_sum_approx": float(sum(cs)),
                      "accepted_state_pair_acoustic_dt_approx": pair_dt,
                      "pair_dt_over_stored_K4_dt": pair_dt/float(dt[i])})
    return {"initial": str(initial), "snapshot": str(snapshot),
            "snapshot_sha256": digest(snapshot), "initial_sha256": digest(initial),
            "time": meta["t"], "cells": len(dt), "matched_native_log": row,
            "Courant": courant, "gamma": gamma,
            "stored_K4_dt_min": base,
            "stored_K4_dt_quantiles": np.quantile(dt, [0, .01, .1, .5, .9, 1]).tolist(),
            "cells_below_dt_multiples": {str(k): int(np.count_nonzero(dt < k*base))
                                         for k in (2, 4, 16, 64)},
            "disk_cells": int(disk.sum()), "disk_stored_K4_dt_min": float(dt[disk].min()),
            "disk_min_over_global_min": float(dt[disk].min()/base),
            "frozen_ideal_cell_update_ratio": float(len(dt)/np.sum(base/dt)),
            "frozen_power2_cell_update_ratio": float(len(dt)/np.sum(np.exp2(-levels.astype(float)))),
            "bins_by_initial_radius": bins, "limiting_cells": pairs,
            "limiting_pair_rows_truncated": len(selected) < len(limiting),
            "checks": {"fixed_mass_W0_finite_positive": True,
                       "complete_unique_IDs_and_native_layout": True,
                       "matching_native_step_time_dt_and_final_GPU_min": True},
            "scope": "Stored last-K4 hydro/geometric dt, not recomputed accepted-state RHS. "
                     "Pair formula uses accepted positions/velocities/EOS and is approximate; "
                     "it does not establish the full winning limiter or native face CSR. "
                     "Frozen cell-update ratios exclude mesh rebuild, communication, "
                     "conservative asynchronous fluxes, limiter evolution and extra global "
                     "bounds; they are not measured or promised runtime speedups."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("initial", type=Path)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("log", type=Path)
    parser.add_argument("--step", type=int, required=True)
    parser.add_argument("--courant", type=float, required=True)
    parser.add_argument("--gamma", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not (np.isfinite(args.courant) and args.courant > 0
            and np.isfinite(args.gamma) and args.gamma > 1 and args.step > 0):
        parser.error("Expected positive finite Courant, gamma>1, step>0")
    result = analyze(args.initial, args.snapshot, args.log, args.step, args.courant, args.gamma)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))

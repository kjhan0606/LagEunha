"""Read-only periodic-seam diagnostic for native Kepler snapshots."""

import argparse
import json
from pathlib import Path
import struct
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot

GAMMA_MINUS_ONE = float(np.float32(1.666667))-1.


def ordered(path):
    with path.open("rb") as stream:
        magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
    if (magic, stride, stress) != (-1, 576, 1) or path.stat().st_size != 16+n*stride+52:
        raise ValueError(f"Incomplete or incompatible snapshot {path}")
    p, meta = read_snapshot(path)
    ids = (p["u4if"].astype(np.uint64) & np.uint64((1 << 56)-1)).astype(np.int64)
    if not np.array_equal(np.sort(ids), np.arange(n)):
        raise ValueError(f"Invalid cell IDs {path}")
    return p[np.argsort(ids)], meta


def weighted_mean(values, weights):
    return float(np.sum(values*weights, dtype=np.longdouble)/
                 np.sum(weights, dtype=np.longdouble))


def stats(mask, density_change, pressure_change, entropy, mass):
    if not np.any(mask):
        return None
    w = mass[mask]
    return {
        "cells": int(np.sum(mask)),
        "density_change_abs_mass_weighted": weighted_mean(np.abs(density_change[mask]), w),
        "density_change_abs_p95": float(np.quantile(np.abs(density_change[mask]), .95)),
        "pressure_change_abs_mass_weighted": weighted_mean(np.abs(pressure_change[mask]), w),
        "pressure_change_abs_p95": float(np.quantile(np.abs(pressure_change[mask]), .95)),
        "entropy_log_signed_mass_weighted": weighted_mean(entropy[mask], w),
        "entropy_log_abs_mass_weighted": weighted_mean(np.abs(entropy[mask]), w),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("endpoint", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    p0, m0 = ordered(args.initial)
    p, mt = ordered(args.endpoint)
    if len(p0) != len(p) or float(m0["t"]) != 0 or float(mt["t"]) <= 0:
        raise ValueError("Expected matching t=0 and later native snapshots")
    lx, ly = float(mt["Lx"]), float(mt["Ly"])
    if (lx, ly) != (6., 6.):
        raise ValueError("This audit expects the six-by-six Kepler box")
    if not np.array_equal(p0["mass"], p["mass"]):
        raise ValueError("Cell masses changed")
    if not (np.all(p0["w2"] == 0) and np.all(p["w2"] == 0)):
        raise ValueError("Expected W=0")
    if not np.all((p["volume"] > 0) & (p["ie"] > 0) & (p["den"] > 0) &
                  (p0["volume"] > 0) & (p0["ie"] > 0) & (p0["den"] > 0)):
        raise ValueError("Nonpositive state")
    xy = np.column_stack((p["x"], p["y"]))
    if not np.all((xy >= 0) & (xy < 6)):
        raise ValueError("Current generator outside periodic box")
    r0 = np.hypot(p0["x"]-3., p0["y"]-3.)
    rt = np.hypot(p["x"]-3., p["y"]-3.)
    edge = np.minimum.reduce((p["x"], p["y"], 6.-p["x"], 6.-p["y"]))
    mass = p0["mass"].astype(float)
    density_change = p["den"].astype(float)/p0["den"].astype(float)-1.
    pressure_change = p["pressure"].astype(float)/p0["pressure"].astype(float)-1.
    entropy = (np.log(p["ie"].astype(float)/p0["ie"].astype(float)) +
               GAMMA_MINUS_ONE*np.log(p["volume"].astype(float)/
                                      p0["volume"].astype(float)))
    if not all(np.isfinite(q).all() for q in (density_change, pressure_change, entropy)):
        raise ValueError("Nonfinite relative field")
    volume = p["volume"].astype(float)
    tree = cKDTree(xy, boxsize=6.)
    candidates = tree.query_pairs(.2*np.sqrt(np.max(volume)), output_type="ndarray")
    if len(candidates):
        ia, ib = candidates[:, 0], candidates[:, 1]
        raw = xy[ia]-xy[ib]
        displacement = (raw+3.)%6.-3.
        d = np.linalg.norm(displacement, axis=1)
        cutoff = .2*np.minimum(np.sqrt(volume[ia]), np.sqrt(volume[ib]))
        short = d < cutoff
        pairs = candidates[short]
        pair_raw = raw[short]
        pair_d = d[short]
    else:
        pairs = np.empty((0, 2), dtype=np.intp)
        pair_raw = np.empty((0, 2))
        pair_d = np.empty(0)
    endpoints = np.unique(pairs)
    velocity = np.column_stack((p["vx"], p["vy"]))
    dv = np.linalg.norm(velocity[pairs[:, 0]]-velocity[pairs[:, 1]], axis=1)
    result = {
        "initial": str(args.initial), "endpoint": str(args.endpoint),
        "time": float(mt["t"]), "cells": len(p),
        "short_pairs": len(pairs), "short_pair_endpoint_cells": len(endpoints),
        "short_pairs_crossing_periodic_seam":
            int(np.sum(np.any(np.abs(pair_raw) > 3., axis=1))),
        "short_pair_distance_minmax":
            [float(np.min(pair_d)), float(np.max(pair_d))] if len(pairs) else None,
        "short_pair_relative_speed_min_median_max":
            [float(np.min(dv)), float(np.median(dv)), float(np.max(dv))]
            if len(pairs) else None,
        "short_pair_current_radius_minmax":
            [float(np.min(rt[endpoints])), float(np.max(rt[endpoints]))]
            if len(endpoints) else None,
        "short_pair_edge_distance_min_median_max":
            [float(np.min(edge[endpoints])), float(np.median(edge[endpoints])),
             float(np.max(edge[endpoints]))] if len(endpoints) else None,
        "regions": {
            "disk_initial_0p5_to_2": stats((r0 >= .5)&(r0 < 2.), density_change,
                                            pressure_change, entropy, mass),
            "outer_initial_r_ge_3": stats(r0 >= 3., density_change,
                                           pressure_change, entropy, mass),
            "current_edge_lt_0p2": stats(edge < .2, density_change,
                                          pressure_change, entropy, mass),
            "current_edge_gt_0p5": stats(edge > .5, density_change,
                                          pressure_change, entropy, mass),
        },
        "checks": {
            "mass_fixed_exact": bool(np.array_equal(p0["mass"], p["mass"])),
            "W0_exact": bool(np.all(p0["w2"] == 0) and np.all(p["w2"] == 0)),
            "density_mass_volume_max_relative_below_1e_minus_6": bool(
                np.max(np.abs(p["den"].astype(float)/(mass/volume)-1.)) < 1e-6),
        },
        "scope": "Native snapshot periodic generator-distance and radial/edge field audit. Short pairs are below the restoring-pressure distance cutoff, but pressure-work activation is determined separately by the native force loop. Region statistics are observational, not a causality proof.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))
    if not all(result["checks"].values()):
        raise SystemExit("Seam checkpoint gate failed")


if __name__ == "__main__":
    main()

"""Periodic short-pair screen for the polar Kepler layout.

This screen uses generator distances and the pair-pressure cutoff.
It is not a face-level pressure audit. Exact-orbit predictions for initial
r >= half the box width are not physical orbital controls after wrapping.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from audit_kepler_polar_checkpoint_20261026 import CX, CY, EPS, load_ordered

BINS = ((0., .1), (.1, .5), (.5, 2.), (2., 2.35), (2.35, 3.), (3., 5.))


def screen(xy, volumes, r0, box_size):
    tree = cKDTree(xy, boxsize=box_size)
    distance, nearest = tree.query(xy, k=2, workers=1)
    distance, nearest = distance[:, 1], nearest[:, 1]
    cutoff = .2 * np.minimum(np.sqrt(volumes), np.sqrt(volumes[nearest]))
    if not np.all(np.isfinite(cutoff) & (cutoff > 0)):
        raise ValueError("Invalid pair-pressure cutoff")
    ratio = distance / cutoff
    near = ratio < 1.
    largest_cutoff = .2 * np.sqrt(np.max(volumes))
    if largest_cutoff >= box_size / 4:
        raise ValueError("Unphysically large pair cutoff for a periodic screen")
    candidate = tree.query_pairs(largest_cutoff, output_type="ndarray")
    if len(candidate):
        ia, ib = candidate[:, 0], candidate[:, 1]
        delta = (xy[ia] - xy[ib] + box_size/2) % box_size - box_size/2
        pair_distance = np.linalg.norm(delta, axis=1)
        pair_cutoff = .2 * np.minimum(np.sqrt(volumes[ia]), np.sqrt(volumes[ib]))
        active_distance = pair_distance < pair_cutoff
        short_pair_rmin = np.minimum(r0[ia[active_distance]], r0[ib[active_distance]])
        short_endpoints = np.unique(candidate[active_distance])
    else:
        active_distance = np.empty(0, dtype=bool)
        short_pair_rmin = np.empty(0)
        short_endpoints = np.empty(0, dtype=np.intp)
    current_radius = np.hypot(xy[:, 0]-CX, xy[:, 1]-CY)
    return {
        "near_cells": int(np.sum(near)),
        "short_generator_pairs": int(np.sum(active_distance)),
        "short_pair_endpoint_cells": int(len(short_endpoints)),
        "short_pair_endpoint_initial_radius_minmax":
            [float(np.min(r0[short_endpoints])), float(np.max(r0[short_endpoints]))]
            if len(short_endpoints) else None,
        "short_pair_endpoint_current_radius_minmax":
            [float(np.min(current_radius[short_endpoints])),
             float(np.max(current_radius[short_endpoints]))]
            if len(short_endpoints) else None,
        "minimum_distance_over_cutoff": float(np.min(ratio)),
        "radial_bins": [
            {"rlo": lo, "rhi": hi, "cells": int(np.sum(sel)),
             "near_cells": int(np.sum(near & sel)),
             "short_generator_pairs": int(np.sum((short_pair_rmin >= lo) &
                                                   (short_pair_rmin < hi))),
             "minimum_distance_over_cutoff": float(np.min(ratio[sel]))}
            for lo, hi in BINS
            if np.any(sel := ((r0 >= lo) & (r0 < hi)))
        ],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("--endpoint", type=Path)
    parser.add_argument("--times", nargs="*", type=float, default=[])
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    p0, m0 = load_ordered(args.initial)
    if float(m0["t"]) != 0 or float(m0["Lx"]) != float(m0["Ly"]):
        raise ValueError("Expected square-domain initial snapshot at t=0")
    box_size = float(m0["Lx"])
    x0, y0 = p0["x"].astype(float)-CX, p0["y"].astype(float)-CY
    r0 = np.hypot(x0, y0)
    volumes0 = p0["volume"].astype(float)
    outputs = []
    for t in args.times:
        phase = t * (r0*r0 + EPS*EPS)**(-.75)
        xy = np.column_stack((np.mod(CX+x0*np.cos(phase)-y0*np.sin(phase), box_size),
                              np.mod(CY+x0*np.sin(phase)+y0*np.cos(phase), box_size)))
        outputs.append(dict(kind="ballistic_with_initial_cell_volumes", time=t,
                            **screen(xy, volumes0, r0, box_size)))
    if args.endpoint is not None:
        pt, mt = load_ordered(args.endpoint)
        xy = np.column_stack((pt["x"], pt["y"]))
        if not np.array_equal(pt["mass"], p0["mass"]):
            raise ValueError("Native checkpoint changed cell masses")
        outputs.append(dict(kind="native_with_native_cell_volumes", time=float(mt["t"]),
                            **screen(xy, pt["volume"].astype(float), r0, box_size)))
    result = {
        "initial": str(args.initial), "endpoint": str(args.endpoint) if args.endpoint else None,
        "distance_cutoff": "d_c=0.2 min(sqrt(V_i),sqrt(V_j))",
        "periodic_box_size": box_size,
        "screens": outputs,
        "scope": "Periodic generator-distance screen, not a Voronoi-face count or pressure-work measurement. Ballistic entries retain initial volumes as a proxy. Circular orbits crossing the periodic boundary at initial r>=3 are not native dynamical solutions.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

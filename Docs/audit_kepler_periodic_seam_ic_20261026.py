"""Measure the initial Kepler velocity jump across the periodic box seams.

This is an IC diagnostic, not an evolution or a flux-work calculation.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from audit_kepler_seam_checkpoint_20261026 import ordered


def summarize(values):
    return {
        "min": float(np.min(values)),
        "median": float(np.median(values)),
        "p95": float(np.quantile(values, 0.95)),
        "max": float(np.max(values)),
    }


def seam_pairs(particles, axis, box_size=6.0, strip_width=0.04, cutoff=0.026):
    xy = np.column_stack((particles["x"], particles["y"])).astype(float)
    velocity = np.column_stack((particles["vx"], particles["vy"])).astype(float)
    low = np.flatnonzero(xy[:, axis] < strip_width)
    high = np.flatnonzero(xy[:, axis] > box_size - strip_width)
    shifted_high = xy[high].copy()
    shifted_high[:, axis] -= box_size
    tree = cKDTree(shifted_high)
    neighbours = tree.query_ball_point(xy[low], cutoff)
    pairs = [(int(i), int(high[j])) for i, near in zip(low, neighbours) for j in near]
    if not pairs:
        raise ValueError(f"No pairs across axis {axis} seam")
    pair = np.asarray(pairs, dtype=np.int64)
    i, j = pair[:, 0], pair[:, 1]
    displacement = xy[i] - xy[j]
    displacement[:, axis] += box_size
    distance = np.linalg.norm(displacement, axis=1)
    dv = velocity[i] - velocity[j]
    tangent = 1 - axis
    gamma = 5.0 / 3.0
    cs_i = np.sqrt(gamma * particles["pressure"][i] / particles["den"][i])
    cs_j = np.sqrt(gamma * particles["pressure"][j] / particles["den"][j])
    pair_cs = 0.5 * (cs_i + cs_j)
    return {
        "pairs_with_distance_below_cutoff": len(pair),
        "distance": summarize(distance),
        "absolute_tangential_velocity_jump": summarize(np.abs(dv[:, tangent])),
        "absolute_normal_velocity_jump": summarize(np.abs(dv[:, axis])),
        "absolute_mean_normal_speed": summarize(np.abs(0.5 * (velocity[i, axis] + velocity[j, axis]))),
        "tangential_jump_over_local_sound_speed": summarize(np.abs(dv[:, tangent]) / pair_cs),
        "sound_speed": summarize(pair_cs),
        "strip_width": strip_width,
        "pair_distance_cutoff": cutoff,
    }


def analytic_jump(offset, eps=0.01175):
    r2 = 9.0 + offset * offset + eps * eps
    radius = np.sqrt(9.0 + offset * offset)
    rho = 0.01 + (1.0 + (radius - 2.0) / 0.1) ** -3
    cs = np.sqrt((5.0 / 3.0) * 1e-6 / rho)
    velocity_jump = 6.0 / r2 ** 0.75
    normal_acceleration_jump = 6.0 / r2 ** 1.5
    return {
        "offset_from_midline": offset,
        "velocity_tangent_jump": velocity_jump,
        "normal_acceleration_jump": normal_acceleration_jump,
        "local_sound_speed": cs,
        "tangential_jump_over_sound_speed": velocity_jump / cs,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("snapshots", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = {
        "initial_snapshots": [],
        "analytic_centered_seam": [analytic_jump(y) for y in (0.0, 1.0, 2.0, 3.0)],
        "interpretation": (
            "The analytic centered Kepler velocity and acceleration are not smooth on the "
            "periodically identified 6x6 seams. A tangential jump with nonzero normal "
            "advection away from the seam midline is not a steady circular-flow IC. "
            "The finite-distance native pairs measure this mismatch but do not prove "
            "its later effect on the disk."
        ),
    }
    for path in args.snapshots:
        p, meta = ordered(path)
        if float(meta["t"]) != 0.0 or (float(meta["Lx"]), float(meta["Ly"])) != (6.0, 6.0):
            raise ValueError(f"Expected a t=0 snapshot in a 6x6 box: {path}")
        if not (np.all(p["w2"] == 0) and np.all(p["den"] > 0) and
                np.all(p["pressure"] > 0)):
            raise ValueError(f"Unexpected initial state: {path}")
        result["initial_snapshots"].append({
            "path": str(path),
            "cells": len(p),
            "x_seam": seam_pairs(p, 0),
            "y_seam": seam_pairs(p, 1),
        })
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

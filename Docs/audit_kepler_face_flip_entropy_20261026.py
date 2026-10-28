"""Correlate Kepler disk Voronoi face turnover with native cell entropy error."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi
from scipy.stats import spearmanr

from audit_kepler_polar_checkpoint_20261026 import (
    CX, CY, GAMMA_MINUS_ONE, load_ordered,
)


def disk_face_keys(p, disk):
    xy = np.column_stack((p["x"], p["y"]))
    vor = Voronoi(xy, qhull_options="Qbb Qc Qz")
    pairs = vor.ridge_points
    pairs = pairs[disk[pairs[:, 0]] & disk[pairs[:, 1]]]
    n = len(p)
    a = np.minimum(pairs[:, 0], pairs[:, 1]).astype(np.uint64)
    b = np.maximum(pairs[:, 0], pairs[:, 1]).astype(np.uint64)
    keys = np.sort(a*np.uint64(n)+b)
    if len(np.unique(keys)) != len(keys):
        raise ValueError("Duplicate disk face")
    return keys


def weighted_mean(values, weights):
    return float(np.sum(values*weights, dtype=np.longdouble) /
                 np.sum(weights, dtype=np.longdouble))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("endpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    p0, m0 = load_ordered(args.initial)
    p, mt = load_ordered(args.endpoint)
    if float(m0["t"]) != 0 or float(mt["t"]) <= 0:
        raise ValueError("Expected t=0 and a later complete native checkpoint")
    if not np.array_equal(p0["mass"], p["mass"]):
        raise ValueError("Cell masses changed")
    if not (np.all(p0["w2"] == 0) and np.all(p["w2"] == 0)):
        raise ValueError("Expected W=0 Voronoi geometry")
    r0 = np.hypot(p0["x"]-CX, p0["y"]-CY)
    disk = (r0 >= .5) & (r0 < 2.)
    initial_keys = disk_face_keys(p0, disk)
    final_keys = disk_face_keys(p, disk)
    lost = np.setdiff1d(initial_keys, final_keys, assume_unique=True)
    added = np.setdiff1d(final_keys, initial_keys, assume_unique=True)
    n = len(p)
    changes = np.zeros(n, dtype=np.int32)
    for keys in (lost, added):
        first = (keys // np.uint64(n)).astype(np.intp)
        second = (keys % np.uint64(n)).astype(np.intp)
        changes += np.bincount(first, minlength=n).astype(np.int32)
        changes += np.bincount(second, minlength=n).astype(np.int32)

    entropy = (np.log(p["ie"].astype(float)/p0["ie"].astype(float)) +
               GAMMA_MINUS_ONE*np.log(p["volume"].astype(float)/
                                      p0["volume"].astype(float)))
    density = p["den"].astype(float)/p0["den"].astype(float)-1.
    if not np.all(np.isfinite(entropy[disk]) & np.isfinite(density[disk])):
        raise ValueError("Invalid disk fields")
    mass = p0["mass"].astype(float)
    groups = []
    for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 5), (5, 10000)):
        selected = disk & (changes >= lo) & (changes < hi)
        if not np.any(selected):
            continue
        groups.append({
            "changed_face_ends_min": lo, "changed_face_ends_max_exclusive": hi,
            "cells": int(np.sum(selected)),
            "entropy_log_abs_mass_weighted": weighted_mean(np.abs(entropy[selected]), mass[selected]),
            "entropy_log_signed_mass_weighted": weighted_mean(entropy[selected], mass[selected]),
            "density_change_abs_mass_weighted": weighted_mean(np.abs(density[selected]), mass[selected]),
        })
    radial = []
    for lo, hi in ((.5, .75), (.75, 1.), (1., 1.5), (1.5, 2.)):
        selected = disk & (r0 >= lo) & (r0 < hi)
        changed = selected & (changes > 0)
        unchanged = selected & (changes == 0)
        if not (np.any(changed) and np.any(unchanged)):
            continue
        radial.append({
            "rlo": lo, "rhi": hi, "cells": int(np.sum(selected)),
            "changed_cells": int(np.sum(changed)),
            "changed_face_ends_vs_abs_entropy_spearman": float(
                spearmanr(changes[selected], np.abs(entropy[selected])).statistic),
            "changed_face_ends_vs_abs_density_change_spearman": float(
                spearmanr(changes[selected], np.abs(density[selected])).statistic),
            "abs_entropy_changed_mass_weighted": weighted_mean(
                np.abs(entropy[changed]), mass[changed]),
            "abs_entropy_unchanged_mass_weighted": weighted_mean(
                np.abs(entropy[unchanged]), mass[unchanged]),
        })
    result = {
        "initial": str(args.initial), "endpoint": str(args.endpoint),
        "time": float(mt["t"]), "disk_cells": int(np.sum(disk)),
        "initial_disk_faces": int(len(initial_keys)),
        "final_disk_faces": int(len(final_keys)),
        "lost_faces": int(len(lost)), "new_faces": int(len(added)),
        "cells_with_changed_disk_faces": int(np.sum(disk & (changes > 0))),
        "changed_face_ends_vs_abs_entropy_spearman": float(
            spearmanr(changes[disk], np.abs(entropy[disk])).statistic),
        "changed_face_ends_vs_abs_density_change_spearman": float(
            spearmanr(changes[disk], np.abs(density[disk])).statistic),
        "groups": groups,
        "radial_bins": radial,
        "scope": "Disk-disk Voronoi face turnover versus per-cell native entropy/density error; association is not causal face-flux attribution. Disk-boundary faces with atmosphere cells are excluded.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

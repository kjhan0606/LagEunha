"""Separate native sparse-disk density change from exact-orbit Voronoi geometry.

This is a read-only diagnostic of an existing accepted snapshot, not a gas
simulation or a substitute for the energy/entropy/positivity gate.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import IC, digest, snapshot_time


EPS = .01175
BINS = ((.5, .75), (.75, 1.), (1., 2.), (.5, 2.))


def weighted_abs(values, weights):
    return float(np.sum(np.abs(values) * weights, dtype=np.longdouble) /
                 np.sum(weights, dtype=np.longdouble))


def disk_geometry(p0, t, selected):
    x0 = p0["x"].astype(float) - 3.
    y0 = p0["y"].astype(float) - 3.
    r0 = np.hypot(x0, y0)
    phase = t * (r0 * r0 + EPS * EPS) ** (-.75)
    xy = np.column_stack((3. + x0 * np.cos(phase) - y0 * np.sin(phase),
                          3. + x0 * np.sin(phase) + y0 * np.cos(phase)))
    vor = Voronoi(xy, qhull_options="Qbb Qc Qz")
    area = np.empty(len(selected), dtype=float)
    radius = np.empty(len(selected), dtype=float)
    for k, i in enumerate(selected):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded exact-orbit disk cell {i}")
        vertices = vor.vertices[region]
        area[k] = .5 * abs(np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1)) -
                             np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1)))
        radius[k] = np.max(np.linalg.norm(vertices - xy[i], axis=1))
    edge = np.min(np.column_stack((xy[selected, 0], xy[selected, 1],
                                   6. - xy[selected, 0], 6. - xy[selected, 1])),
                  axis=1)
    if not np.all(np.isfinite(area) & (area > 0) &
                  np.isfinite(radius) & (radius > 0) & (edge > 2 * radius)):
        raise ValueError("Exact-orbit disk polygon fails periodic-image safety")
    return xy, area, float(np.min(edge / (2 * radius)))


def analyze(initial, endpoint, count, target):
    if not np.isfinite(target) or target <= 0:
        raise ValueError("Target time must be finite and positive")
    if abs(snapshot_time(initial, count)) > 1e-12 or abs(
            snapshot_time(endpoint, count) - target) > 1e-10:
        raise ValueError("Incomplete or wrong-time native checkpoint")
    p0, m0 = ordered(initial, count)
    p, m = ordered(endpoint, count)
    if abs(float(m0["t"])) > 1e-12 or abs(float(m["t"]) - target) > 1e-10:
        raise ValueError("Header/footer time mismatch")
    fields = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2")
    finite = bool(all(np.isfinite(q[key]).all() for q in (p0, p) for key in fields))
    if not finite:
        raise ValueError("Nonfinite native state")
    initial_area = float(np.sum(p0["volume"].astype(float), dtype=np.longdouble))
    final_area = float(np.sum(p["volume"].astype(float), dtype=np.longdouble))
    rho_error = max(float(np.max(np.abs(
        q["den"].astype(float) /
        (q["mass"].astype(float) / q["volume"].astype(float)) - 1.)))
        for q in (p0, p))
    checks = {
        "finite": finite,
        "fixed_mass_exact": bool(np.array_equal(p0["mass"], p["mass"])),
        "W0_exact": bool(np.all(p0["w2"] == 0) and np.all(p["w2"] == 0)),
        "positive_mass_volume_density": bool(all(np.all(q[key] > 0)
                                                  for q in (p0, p)
                                                  for key in ("mass", "volume", "den"))),
        "periodic_area_tiles_box": bool(abs(initial_area - 36.) < 1e-6 and
                                        abs(final_area - 36.) < 1e-5),
        "rho_equals_mass_over_volume": bool(np.isfinite(rho_error) and rho_error < 1e-6),
    }
    if not all(checks.values()):
        raise ValueError(f"Native state fails geometry prerequisites: {checks}")
    x0 = p0["x"].astype(float) - 3.
    y0 = p0["y"].astype(float) - 3.
    r0 = np.hypot(x0, y0)
    selected = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    if not len(selected):
        raise ValueError("Missing initial Kepler disk")
    xy, ballistic_area, safety = disk_geometry(p0, target, selected)
    mass = p0["mass"][selected].astype(float)
    initial_volume = p0["volume"][selected].astype(float)
    final_volume = p["volume"][selected].astype(float)
    initial_density = p0["den"][selected].astype(float)
    native_density = p["den"][selected].astype(float)
    native_change = native_density / initial_density - 1.
    ballistic_change = initial_volume / ballistic_area - 1.
    # This compares densities, not the difference between two absolute norms.
    native_vs_ballistic = (mass / final_volume) / (mass / ballistic_area) - 1.
    dx = (p["x"][selected].astype(float) - xy[selected, 0] + 3.) % 6. - 3.
    dy = (p["y"][selected].astype(float) - xy[selected, 1] + 3.) % 6. - 3.
    position_error = np.hypot(dx, dy)
    rows = []
    for lo, hi in BINS:
        take = (r0[selected] >= lo) & (r0[selected] < hi)
        w = mass[take]
        rows.append({
            "rlo": lo, "rhi": hi, "cells": int(np.sum(take)),
            "native_density_change_abs_mass_weighted": weighted_abs(native_change[take], w),
            "ballistic_density_change_abs_mass_weighted": weighted_abs(ballistic_change[take], w),
            "native_minus_ballistic_density_abs_mass_weighted": weighted_abs(native_vs_ballistic[take], w),
            "native_minus_ballistic_density_p95": float(np.quantile(np.abs(native_vs_ballistic[take]), .95)),
            "native_minus_ballistic_density_max": float(np.max(np.abs(native_vs_ballistic[take]))),
            "position_error_abs_mass_weighted": weighted_abs(position_error[take], w),
            "position_error_p95": float(np.quantile(position_error[take], .95)),
        })
    return {
        "initial": str(initial.resolve()), "initial_sha256": digest(initial),
        "endpoint": str(endpoint.resolve()), "endpoint_sha256": digest(endpoint),
        "time": target, "n": count, "disk_cells": len(selected),
        "initial_area_sum": initial_area, "endpoint_area_sum": final_area,
        "rho_mass_volume_max_relative": rho_error,
        "minimum_periodic_image_safety_ratio": safety, "checks": checks,
        "bins": rows,
        "scope": "Native fixed-mass W0 density versus exact softened circular-orbit Voronoi geometry; does not audit entropy, energy, heating, or long-time stability",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("endpoint", type=Path)
    parser.add_argument("--expected-count", type=int, choices=tuple(IC), required=True)
    parser.add_argument("--time", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = analyze(args.initial, args.endpoint, args.expected_count, args.time)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"time": result["time"], "n": result["n"],
                      "disk": result["bins"][-1]}, indent=2))


if __name__ == "__main__":
    main()

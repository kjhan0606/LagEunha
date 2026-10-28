"""Audit any complete native polar Kepler checkpoint against exact-orbit geometry.

Read-only inputs. The output is a new JSON receipt; the active run is never
modified. Disk cells are far enough from the 6x6 edge that omitted periodic
images cannot affect their Voronoi polygons at these early checkpoints.
"""

import argparse
import json
import struct
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi

SNAPSHOT_MODULE = Path("/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
sys.path.insert(0, str(SNAPSHOT_MODULE))
from audit_snapshot import read_snapshot

N = 384 * 384
EPS = .01175
CX = CY = 3.0
GAMMA_MINUS_ONE = float(np.float32(1.666667)) - 1.0


def load_ordered(path: Path):
    with path.open("rb") as stream:
        magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
    if (magic, n, stride, stress) != (-1, N, 576, 1):
        raise ValueError(f"Unexpected native snapshot header: {path}")
    if path.stat().st_size != 16 + n * stride + 52:
        raise ValueError(f"Incomplete or overlong native snapshot: {path}")
    data, meta = read_snapshot(path)
    ids = (data["u4if"].astype(np.uint64) &
           np.uint64((1 << 56) - 1)).astype(np.int64)
    if not np.array_equal(np.sort(ids), np.arange(N)):
        raise ValueError(f"Snapshot cell IDs incomplete: {path}")
    return data[np.argsort(ids)], meta


def weighted_mean(values, weights):
    return float(np.sum(values * weights, dtype=np.longdouble) /
                 np.sum(weights, dtype=np.longdouble))


def disk_ballistic_areas(p0, t, selected):
    x0, y0 = p0["x"] - CX, p0["y"] - CY
    r0 = np.hypot(x0, y0)
    angle = t * (r0*r0 + EPS*EPS)**(-.75)
    xy = np.column_stack((CX + x0*np.cos(angle) - y0*np.sin(angle),
                          CY + x0*np.sin(angle) + y0*np.cos(angle)))
    vor = Voronoi(xy, qhull_options="Qbb Qc Qz")
    areas = np.empty(len(selected), dtype=np.float64)
    for k, i in enumerate(selected):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded disk cell in exact-orbit control: {i}")
        vertices = vor.vertices[region]
        areas[k] = .5 * abs(np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1)) -
                               np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1)))
    if not np.all(np.isfinite(areas) & (areas > 0)):
        raise ValueError("Invalid ballistic disk area")
    return xy, areas


def analyze(initial: Path, endpoint: Path):
    p0, meta0 = load_ordered(initial)
    p, meta = load_ordered(endpoint)
    t = float(meta["t"])
    if meta0["t"] != 0 or not (t > 0 and t <= .825 + 1e-12):
        raise ValueError("Expected exact t=0 initial and early endpoint")
    finite = all(np.isfinite(q[field]).all() for q in (p0, p) for field in
                 ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2"))
    checks = {"finite": bool(finite),
              "mass_exact": bool(np.array_equal(p0["mass"], p["mass"])),
              "W0_exact": bool(np.all(p0["w2"] == 0) and np.all(p["w2"] == 0)),
              "volume_positive": bool(np.all(p0["volume"] > 0) and np.all(p["volume"] > 0)),
              "internal_energy_positive": bool(np.all(p0["ie"] > 0) and np.all(p["ie"] > 0))}
    if not all(checks.values()):
        raise ValueError(f"Native state failed physical gate: {checks}")

    x0, y0 = p0["x"] - CX, p0["y"] - CY
    r0 = np.hypot(x0, y0)
    disk = (r0 >= .5) & (r0 < 2.)
    selected = np.flatnonzero(disk)
    xy_ballistic, area_ballistic = disk_ballistic_areas(p0, t, selected)
    position_error = np.hypot(p["x"] - xy_ballistic[:, 0],
                              p["y"] - xy_ballistic[:, 1])
    native_density_change = p["den"].astype(float) / p0["den"].astype(float) - 1
    pressure_change = p["pressure"].astype(float) / p0["pressure"].astype(float) - 1
    entropy_ratio = p["ie"].astype(float) / (
        p0["ie"].astype(float) *
        (p0["volume"].astype(float) / p["volume"].astype(float))**GAMMA_MINUS_ONE)
    geom_density_change = p0["volume"][selected].astype(float) / area_ballistic - 1
    native_density = p["mass"][selected].astype(float) / p["volume"][selected].astype(float)
    ballistic_density = p0["mass"][selected].astype(float) / area_ballistic
    geom_residual = native_density / ballistic_density - 1

    weights = p0["mass"].astype(float)
    bins = []
    for lo, hi in ((0., .1), (.1, .5), (.5, .75), (.75, 1.), (1., 2.), (.5, 2.), (2., 2.35)):
        mask = (r0 >= lo) & (r0 < hi)
        if not np.any(mask):
            continue
        local = (r0[selected] >= lo) & (r0[selected] < hi)
        row = {"rlo": lo, "rhi": hi, "cells": int(np.sum(mask)),
               "native_density_abs_mass_weighted": weighted_mean(abs(native_density_change[mask]), weights[mask]),
               "pressure_abs_mass_weighted": weighted_mean(abs(pressure_change[mask]), weights[mask]),
               "entropy_abs_mass_weighted": weighted_mean(abs(entropy_ratio[mask] - 1), weights[mask]),
               "position_error_mass_weighted": weighted_mean(position_error[mask], weights[mask]),
               "position_error_max": float(np.max(position_error[mask])),
               "angular_momentum_relative": float(
                   np.sum(weights[mask] * ((p["x"][mask]-CX)*p["vy"][mask] -
                                           (p["y"][mask]-CY)*p["vx"][mask])) /
                   np.sum(weights[mask] * (x0[mask]*p0["vy"][mask] - y0[mask]*p0["vx"][mask])) - 1)}
        if np.any(local):
            w = weights[selected][local]
            row.update(ballistic_geometry_density_abs_mass_weighted=weighted_mean(abs(geom_density_change[local]), w),
                       native_minus_ballistic_density_abs_mass_weighted=weighted_mean(abs(geom_residual[local]), w),
                       native_minus_ballistic_density_p95=float(np.quantile(abs(geom_residual[local]), .95)))
        bins.append(row)
    return {"initial": str(initial), "endpoint": str(endpoint), "time": t,
            "n": N, "checks": checks,
            "minimum_internal_energy": float(np.min(p["ie"])),
            "minimum_volume": float(np.min(p["volume"])),
            "disk_geometry_control": "Independent unweighted Voronoi areas on exact softened circular generator orbits; only initial-r<2 disk is compared",
            "bins": bins,
            "scope": "Native CUDA+MPI+OpenMP W=0 fixed-mass early polar Kepler checkpoint; no long-time stability claim"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("initial", type=Path)
    ap.add_argument("endpoint", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    result = analyze(args.initial, args.endpoint)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

"""Exact circular-orbit Voronoi geometry at late Kepler comparison times.

This is a geometry control on the native initial generators and fixed masses,
not a gas evolution. Only bounded disk polygons certified independent of
omitted periodic images are used.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi

from audit_kepler_polar_checkpoint_20261026 import CX, CY, EPS, load_ordered

TIMES = (0., 6.5, 17.778891990936625, 120., 600., 2221.441469079183)
BINS = ((.5, .75), (.75, 1.), (1., 2.), (.5, 2.))


def mean_abs(values, weights):
    return float(np.sum(np.abs(values)*weights, dtype=np.longdouble) /
                 np.sum(weights, dtype=np.longdouble))


def initial(path):
    p, meta = load_ordered(path)
    if float(meta["t"]) != 0 or not np.all(p["w2"] == 0):
        raise ValueError("Expected complete W=0 initial snapshot")
    x = p["x"].astype(float)-CX
    y = p["y"].astype(float)-CY
    r = np.hypot(x, y)
    selected = np.flatnonzero((r >= .5) & (r < 2.))
    return p, (x, y, r, selected)


def one_time(p, geometry, t):
    x, y, r0, selected = geometry
    phase = t*(r0*r0 + EPS*EPS)**(-.75)
    xy = np.column_stack((CX+x*np.cos(phase)-y*np.sin(phase),
                          CY+x*np.sin(phase)+y*np.cos(phase)))
    vor = Voronoi(xy, qhull_options="Qbb Qc Qz")
    area = np.empty(len(selected), dtype=float)
    circumradius = np.empty(len(selected), dtype=float)
    for k, i in enumerate(selected):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded disk cell at t={t}, i={i}")
        vertices = vor.vertices[region]
        area[k] = .5*abs(np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1))-
                           np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1)))
        circumradius[k] = np.max(np.linalg.norm(vertices-xy[i], axis=1))
    if not np.all(np.isfinite(area) & (area > 0) & np.isfinite(circumradius)):
        raise ValueError(f"Invalid disk Voronoi geometry at t={t}")
    edge = np.min(np.column_stack((xy[selected, 0], xy[selected, 1],
                                   6.-xy[selected, 0], 6.-xy[selected, 1])), axis=1)
    if not np.all(edge > 2*circumradius):
        raise ValueError(f"Disk cell may depend on omitted periodic image at t={t}")
    change = p["volume"][selected].astype(float)/area-1.
    mass = p["mass"][selected].astype(float)
    rows = []
    for lo, hi in BINS:
        local = (r0[selected] >= lo) & (r0[selected] < hi)
        rows.append({"rlo": lo, "rhi": hi, "cells": int(np.sum(local)),
                     "density_change_abs_mass_weighted": mean_abs(change[local], mass[local]),
                     "density_change_abs_p95": float(np.quantile(np.abs(change[local]), .95)),
                     "density_change_abs_max": float(np.max(np.abs(change[local])))})
    return {"time": t, "disk_cells": int(len(selected)),
            "minimum_periodic_image_safety_ratio": float(np.min(edge/(2*circumradius))),
            "bins": rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("q1", type=Path)
    parser.add_argument("q2", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cases = {}
    for name, path in (("q1", args.q1), ("q2", args.q2)):
        p, geometry = initial(path)
        cases[name] = {"initial": str(path), "times": [one_time(p, geometry, t) for t in TIMES]}
    comparison = []
    for q1, q2 in zip(cases["q1"]["times"], cases["q2"]["times"]):
        if q1["time"] != q2["time"]:
            raise ValueError("Mismatched comparison time")
        comparison.append({
            "time": q1["time"],
            "q2_over_q1_disk_density_abs_mass_weighted":
                q2["bins"][-1]["density_change_abs_mass_weighted"] /
                q1["bins"][-1]["density_change_abs_mass_weighted"],
        })
    result = {
        "cases": cases, "comparison": comparison,
        "scope": "Independent exact softened circular-orbit Voronoi geometry of native t=0 W=0 generators. Fixed initial masses, no gas evolution, no prediction of viscous heating or late native stability. Disk polygons pass an omitted-periodic-image safety check.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

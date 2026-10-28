"""Audit narrow/wide native Kepler histories at matched saved times.

All radial masks refer to t=0 generators. The periodic seam is a different
global problem from an isolated Kepler disk; spatial trends are diagnostic,
not a causal proof of boundary contamination or irreversible heating.
"""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


GM1 = float(np.float32(1.666667)) - 1.
EPS = .01175
BINS = ((0., .1), (.1, .2), (.2, .5), (.5, .75), (.75, 1.),
        (1., 2.), (.5, 2.), (2., 2.35), (2.35, 2.75), (2.75, 5.))
STATE_FIELDS = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den",
                "pressure", "w2")


def weighted_mean(values, weights):
    return float(np.sum(values * weights, dtype=np.longdouble) /
                 np.sum(weights, dtype=np.longdouble))


def step(path):
    m = re.fullmatch(r"khout\.(\d+)\.dat", path.name)
    if not m:
        raise ValueError(f"Unexpected native checkpoint name: {path}")
    return int(m.group(1))


def prepare(initial, count):
    p0, meta = ordered(initial, count)
    if (snapshot_time(initial, count) != 0. or float(meta["t"]) != 0. or
            not all(np.isfinite(p0[key]).all() for key in STATE_FIELDS) or
            not np.all(p0["mass"] > 0) or not np.all(p0["volume"] > 0) or
            not np.all(p0["ie"] > 0) or not np.all(p0["w2"] == 0)):
        raise ValueError(f"Invalid complete t=0 native IC: {initial}")
    x0 = p0["x"].astype(float) - 3.
    y0 = p0["y"].astype(float) - 3.
    r0 = np.hypot(x0, y0)
    mass = p0["mass"].astype(float)
    masks = [(r0 >= lo) & (r0 < hi) for lo, hi in BINS]
    disk = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    disk = disk[np.lexsort((p0["y"][disk], p0["x"][disk]))]
    initial_edge = np.minimum.reduce((p0["x"], p0["y"],
                                      6.-p0["x"], 6.-p0["y"])) < .2
    return {"initial": initial, "count": count, "p0": p0, "x0": x0,
            "y0": y0, "r0": r0, "mass": mass, "masks": masks,
            "disk": disk, "initial_edge": initial_edge,
            "sha256": digest(initial)}


def checkpoints(run, count, target, allow_incomplete):
    found, skipped = [], []
    for path in run.glob("khout.*.dat"):
        s = step(path)
        try:
            t = snapshot_time(path, count)
        except ValueError:
            if not allow_incomplete:
                raise
            skipped.append(str(path))
            continue
        if t <= target + 1e-10:
            found.append((t, s, path))
    found.sort()
    if not found or found[0][0] < .825 - 1e-10:
        raise ValueError(f"Missing restart or unexpected early checkpoint: {run}")
    if any(b[0] <= a[0] or b[1] <= a[1] for a, b in zip(found, found[1:])):
        raise ValueError(f"Non-monotone native checkpoint history: {run}")
    if not allow_incomplete and abs(found[-1][0] - target) > 1e-10:
        raise ValueError(f"Missing exact target-time native endpoint: {run}")
    return found, skipped


def region(mask, mass, density, entropy, increment, vr, orbit, j0, j):
    if not np.any(mask):
        return None
    w = mass[mask]
    l0 = np.sum(w*j0[mask], dtype=np.longdouble)
    lt = np.sum(w*j[mask], dtype=np.longdouble)
    return {
        "cells": int(np.sum(mask)),
        "mass": float(np.sum(w, dtype=np.longdouble)),
        "density_abs": weighted_mean(np.abs(density[mask]), w),
        "entropy_signed": weighted_mean(entropy[mask], w),
        "entropy_abs": weighted_mean(np.abs(entropy[mask]), w),
        "entropy_increment_signed": weighted_mean(increment[mask], w),
        "entropy_increment_abs": weighted_mean(np.abs(increment[mask]), w),
        "radial_velocity_abs": weighted_mean(np.abs(vr[mask]), w),
        "orbit_position_error_abs": weighted_mean(orbit[mask], w),
        "angular_momentum_relative": float(lt/l0 - 1.) if abs(l0) > 0 else None,
    }


def analyze(case, item, previous_entropy):
    t, s, path = item
    p, meta = ordered(path, case["count"])
    p0, mass = case["p0"], case["mass"]
    if (abs(float(meta["t"]) - t) > 1e-10 or
            not all(np.isfinite(p[key]).all() for key in STATE_FIELDS) or
            not np.array_equal(p0["mass"], p["mass"]) or
            not np.all(p["w2"] == 0) or not np.all(p["volume"] > 0) or
            not np.all(p["ie"] > 0)):
        raise ValueError(f"Native fixed-mass/W0/positive state gate failed: {path}")
    area = float(np.sum(p["volume"].astype(float), dtype=np.longdouble))
    rho = p["den"].astype(float)
    volume = p["volume"].astype(float)
    u = p["ie"].astype(float)
    rho_error = float(np.max(np.abs(rho/(mass/volume) - 1.)))
    eos_error = float(np.max(np.abs(
        p["pressure"].astype(float)/(GM1*u/volume) - 1.)))
    if abs(area - 36.) >= 1e-5 or rho_error >= 1e-6 or eos_error >= 1e-6:
        raise ValueError(f"Native area/density/EOS gate failed: {path}")
    entropy = (np.log(u/p0["ie"].astype(float)) +
               GM1*np.log(volume/p0["volume"].astype(float)))
    increment = entropy if previous_entropy is None else entropy-previous_entropy
    density = rho/p0["den"].astype(float) - 1.
    x = p["x"].astype(float) - 3.
    y = p["y"].astype(float) - 3.
    r = np.hypot(x, y)
    phase = t*(case["r0"]**2 + EPS**2)**(-.75)
    xb = case["x0"]*np.cos(phase) - case["y0"]*np.sin(phase)
    yb = case["x0"]*np.sin(phase) + case["y0"]*np.cos(phase)
    orbit = np.hypot((x-xb+3.)%6.-3., (y-yb+3.)%6.-3.)
    vx, vy = p["vx"].astype(float), p["vy"].astype(float)
    vr = np.divide(x*vx+y*vy, r, out=np.zeros_like(r), where=r>0)
    j0 = case["x0"]*p0["vy"].astype(float) - case["y0"]*p0["vx"].astype(float)
    j = x*vy-y*vx
    bins = []
    for (lo, hi), mask in zip(BINS, case["masks"]):
        row = region(mask, mass, density, entropy, increment, vr, orbit, j0, j)
        if row is not None:
            row.update(rlo=lo, rhi=hi)
            bins.append(row)
    current_edge = np.minimum.reduce((p["x"], p["y"],
                                      6.-p["x"], 6.-p["y"])) < .2
    inner0, inner = case["r0"] < .5, r < .5
    inward = ~inner0 & inner
    outward = inner0 & ~inner
    row = {
        "checkpoint": str(path.resolve()), "sha256": digest(path),
        "step": s, "time": t, "area_sum": area,
        "rho_mass_volume_max_relative": rho_error,
        "eos_pressure_max_relative": eos_error,
        "minimum_internal_energy": float(np.min(u)),
        "bins": bins,
        "initial_edge_within_0p2": region(case["initial_edge"], mass, density,
                                           entropy, increment, vr, orbit, j0, j),
        "current_edge_within_0p2": region(current_edge, mass, density,
                                           entropy, increment, vr, orbit, j0, j),
        "inner_edge_generator_crossings": {
            "inward_count": int(np.sum(inward)),
            "outward_count": int(np.sum(outward)),
            "net_inward_mass": float(np.sum(mass[inward], dtype=np.longdouble) -
                                     np.sum(mass[outward], dtype=np.longdouble)),
        },
    }
    return row, p, entropy, density


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--narrow-initial", type=Path, required=True)
    ap.add_argument("--narrow-run", type=Path, required=True)
    ap.add_argument("--wide-initial", type=Path, required=True)
    ap.add_argument("--wide-run", type=Path, required=True)
    ap.add_argument("--target-time", type=float, default=10.)
    ap.add_argument("--allow-incomplete", action="store_true")
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if not np.isfinite(args.target_time) or args.target_time <= .825:
        raise ValueError("Expected finite target after the accepted restart")
    cases = {
        "narrow": prepare(args.narrow_initial, 54312),
        "wide": prepare(args.wide_initial, 55931),
    }
    histories = {}
    for name, run in (("narrow", args.narrow_run), ("wide", args.wide_run)):
        histories[name] = checkpoints(run, cases[name]["count"],
                                      args.target_time, args.allow_incomplete)
    a, b = histories["narrow"][0], histories["wide"][0]
    if len(a) != len(b):
        raise ValueError("Narrow/wide saved checkpoint counts differ")
    if any(abs(x[0]-y[0]) > .003 for x, y in zip(a, b)):
        raise ValueError("Narrow/wide saved times are not matched")
    ni, wi = cases["narrow"]["disk"], cases["wide"]["disk"]
    np0, wp0 = cases["narrow"]["p0"], cases["wide"]["p0"]
    if len(ni) != len(wi) or not all(np.array_equal(np0[k][ni], wp0[k][wi])
                                      for k in ("x", "y", "vx", "vy")):
        raise ValueError("Original narrow/wide disk generators differ")
    disk_mass = cases["narrow"]["mass"][ni]
    if np.max(np.abs(cases["wide"]["mass"][wi]/disk_mass-1.)) >= 1e-6:
        raise ValueError("Original narrow/wide disk masses differ")
    rows = {"narrow": [], "wide": []}
    pairs = []
    previous = {"narrow": None, "wide": None}
    for narrow_item, wide_item in zip(a, b):
        nr, npart, ne, nd = analyze(cases["narrow"], narrow_item,
                                   previous["narrow"])
        wr, wpart, we, wd = analyze(cases["wide"], wide_item,
                                   previous["wide"])
        previous = {"narrow": ne, "wide": we}
        rows["narrow"].append(nr)
        rows["wide"].append(wr)
        pairs.append({
            "time_narrow": nr["time"], "time_wide": wr["time"],
            "original_disk_cells": len(ni),
            "matched_density_change_difference_abs": weighted_mean(
                np.abs(nd[ni]-wd[wi]), disk_mass),
            "matched_entropy_log_difference_abs": weighted_mean(
                np.abs(ne[ni]-we[wi]), disk_mass),
            "matched_position_difference_abs": weighted_mean(
                np.hypot(npart["x"][ni].astype(float)-wpart["x"][wi].astype(float),
                         npart["y"][ni].astype(float)-wpart["y"][wi].astype(float)),
                disk_mass),
        })
    result = {
        "target_time": args.target_time,
        "target_reached": abs(a[-1][0]-args.target_time) < 1e-10 and
                          abs(b[-1][0]-args.target_time) < 1e-10,
        "narrow": {"initial": str(args.narrow_initial.resolve()),
                   "initial_sha256": cases["narrow"]["sha256"],
                   "rows": rows["narrow"], "skipped_incomplete": histories["narrow"][1]},
        "wide": {"initial": str(args.wide_initial.resolve()),
                 "initial_sha256": cases["wide"]["sha256"],
                 "rows": rows["wide"], "skipped_incomplete": histories["wide"][1]},
        "matched_original_disk": pairs,
        "scope": "ID-ordered native fixed-mass W0 histories and matched original disk sites. Signed entropy is distinct from irreversible heating; radial/edge trends do not prove propagation. No GIZMO-length claim.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"rows": len(a), "last_time": a[-1][0],
                      "target_reached": result["target_reached"]}))


if __name__ == "__main__":
    main()

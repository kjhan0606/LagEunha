"""Match identical full-q2/sparse-q2 disk generators at native t=.825."""

import argparse
import hashlib
import json
from pathlib import Path
import struct
import sys

import numpy as np

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot


GM1 = float(np.float32(1.666667)) - 1.
BINS = ((.5, .75), (.75, 1.), (1., 2.), (.5, 2.))


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path, n, time):
    with path.open("rb") as stream:
        header = struct.unpack("<4i", stream.read(16))
        stream.seek(-52, 2)
        footer_time = struct.unpack("<d", stream.read(8))[0]
    if (header != (-1, n, 576, 1) or
            path.stat().st_size != 16 + n*576 + 52 or
            abs(footer_time - time) > 1e-12):
        raise ValueError(f"Incomplete or mistimed native state: {path}")
    p, meta = read_snapshot(path)
    ids = (p["u4if"].astype(np.uint64) & np.uint64((1 << 56)-1)).astype(np.int64)
    if (len(p) != n or not np.array_equal(np.sort(ids), np.arange(n)) or
            abs(float(meta["t"]) - time) > 1e-12):
        raise ValueError(f"Incomplete IDs or header/footer time: {path}")
    return p[np.argsort(ids)]


def mean_abs(x, w):
    return float(np.sum(np.abs(x)*w, dtype=np.longdouble) /
                 np.sum(w, dtype=np.longdouble))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("full_initial", type=Path)
    ap.add_argument("full_endpoint", type=Path)
    ap.add_argument("sparse_initial", type=Path)
    ap.add_argument("sparse_endpoint", type=Path)
    ap.add_argument("--sparse-count", type=int, choices=(54312, 55931), default=54312)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    f0 = load(args.full_initial, 384*384, 0.)
    f1 = load(args.full_endpoint, 384*384, .825)
    s0 = load(args.sparse_initial, args.sparse_count, 0.)
    s1 = load(args.sparse_endpoint, args.sparse_count, .825)
    for a, b in ((f0, f1), (s0, s1)):
        if (not np.array_equal(a["mass"], b["mass"]) or
                not np.all(a["w2"] == 0) or not np.all(b["w2"] == 0) or
                not all(np.isfinite(q[k]).all() for q in (a, b)
                        for k in ("x", "y", "vx", "vy", "mass", "volume", "ie", "den")) or
                not np.all(a["volume"] > 0) or not np.all(b["volume"] > 0) or
                not np.all(a["ie"] > 0) or not np.all(b["ie"] > 0)):
            raise ValueError("Native W0/fixed-mass/positive/finite gate failed")
    def selected(p):
        radius = np.hypot(p["x"].astype(float)-3., p["y"].astype(float)-3.)
        ids = np.flatnonzero((radius >= .5) & (radius < 2.))
        return ids[np.lexsort((p["y"][ids], p["x"][ids]))]
    fi, si = selected(f0), selected(s0)
    if len(fi) != 49038 or len(si) != 49038:
        raise ValueError("Full/sparse original disk-site count changed")
    for k in ("x", "y"):
        if not np.array_equal(f0[k][fi], s0[k][si]):
            raise ValueError("Full and sparse disk generators differ initially")
    f0, f1, s0, s1 = f0[fi], f1[fi], s0[si], s1[si]
    w = f0["mass"].astype(float)
    r0 = np.hypot(f0["x"].astype(float)-3., f0["y"].astype(float)-3.)
    def change(p0, p1):
        density = p1["den"].astype(float)/p0["den"].astype(float)-1.
        entropy = (np.log(p1["ie"].astype(float)/p0["ie"].astype(float)) +
                   GM1*np.log(p1["volume"].astype(float)/p0["volume"].astype(float)))
        return density, entropy
    frho, fentropy = change(f0, f1)
    srho, sentropy = change(s0, s1)
    mass_max_relative = float(np.max(np.abs(s0["mass"].astype(float)/w-1.)))
    if mass_max_relative > 1e-6:
        raise ValueError("The matched initial disk masses are not comparable")
    rows = []
    for lo, hi in BINS:
        take = (r0 >= lo) & (r0 < hi)
        mw = w[take]
        rows.append({
            "rlo": lo, "rhi": hi, "cells": int(np.sum(take)),
            "full_density_abs": mean_abs(frho[take], mw),
            "sparse_density_abs": mean_abs(srho[take], mw),
            "matched_density_difference_abs": mean_abs(srho[take]-frho[take], mw),
            "full_entropy_log_abs": mean_abs(fentropy[take], mw),
            "sparse_entropy_log_abs": mean_abs(sentropy[take], mw),
            "matched_entropy_log_difference_abs": mean_abs(sentropy[take]-fentropy[take], mw),
            "full_entropy_log_signed": float(np.sum(fentropy[take]*mw, dtype=np.longdouble)/
                                             np.sum(mw, dtype=np.longdouble)),
            "sparse_entropy_log_signed": float(np.sum(sentropy[take]*mw, dtype=np.longdouble)/
                                               np.sum(mw, dtype=np.longdouble)),
            "fraction_cells_sparse_smaller_abs_entropy": float(np.mean(np.abs(sentropy[take]) <
                                                                 np.abs(fentropy[take]))),
            "endpoint_volume_difference_abs": mean_abs(
                s1["volume"][take].astype(float)/f1["volume"][take].astype(float)-1., mw),
        })
    report = {
        "time": .825, "disk_sites": len(fi), "initial_disk_mass_max_relative": mass_max_relative,
        "initial_disk_velocity_max_absolute": float(max(
            np.max(np.abs(s0[k].astype(float)-f0[k].astype(float))) for k in ("vx", "vy"))),
        "endpoint_generator_position_difference_abs_mass_weighted": mean_abs(
            np.hypot(s1["x"].astype(float)-f1["x"].astype(float),
                     s1["y"].astype(float)-f1["y"].astype(float)), w),
        "bins": rows,
        "inputs_sha256": {str(path.resolve()): digest(path) for path in
                          (args.full_initial, args.full_endpoint,
                           args.sparse_initial, args.sparse_endpoint)},
        "scope": "Matched identical initial disk sites only; does not attribute the difference to timestep, core, or exterior separately",
    }
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(f"matched {len(fi)} disk sites; report {args.output}")


if __name__ == "__main__":
    main()

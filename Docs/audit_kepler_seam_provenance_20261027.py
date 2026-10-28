"""Trace initially exterior Kepler cells toward the disk in native snapshots.

Cell-centre movement is a Lagrangian provenance diagnostic, not a measured
fluid flux or proof of causal propagation through Voronoi faces.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_native_radial_20261026 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


THRESHOLDS = (2.75, 2.35, 2., .5)


def weighted_mean(value, weight):
    return float(np.sum(value * weight, dtype=np.longdouble) /
                 np.sum(weight, dtype=np.longdouble))


def provenance(name, source, r, rt, mass, density_change, entropy_log, vr):
    if not np.any(source):
        raise ValueError(f"Missing {name} provenance group")
    w = mass[source]
    row = {
        "source": name, "cells": int(np.sum(source)),
        "initial_mass": float(np.sum(w, dtype=np.longdouble)),
        "initial_radius_min": float(np.min(r[source])),
        "current_radius_min": float(np.min(rt[source])),
        "current_radius_p05": float(np.quantile(rt[source], .05)),
        "current_radius_median": float(np.median(rt[source])),
        "density_change_abs_mass_weighted": weighted_mean(np.abs(density_change[source]), w),
        "entropy_log_signed_mass_weighted": weighted_mean(entropy_log[source], w),
        "entropy_log_abs_mass_weighted": weighted_mean(np.abs(entropy_log[source]), w),
        "inward_radial_speed_max": float(np.max(np.maximum(-vr[source], 0.))),
        "crossings": [],
    }
    for boundary in THRESHOLDS:
        take = source & (rt < boundary)
        row["crossings"].append({
            "radius": boundary, "cells": int(np.sum(take)),
            "fixed_mass": float(np.sum(mass[take], dtype=np.longdouble)),
        })
    return row


def analyze(initial, endpoint, expected_time):
    p0, m0 = ordered(initial)
    p, m = ordered(endpoint)
    t = float(m["t"])
    if (len(p0) != len(p) or abs(float(m0["t"])) > 1e-12 or
            not np.isfinite(t) or t <= 0 or
            abs(t - expected_time) > 1e-10):
        raise ValueError("Mismatched native snapshot count or time")
    fields = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2")
    finite = bool(all(np.isfinite(q[key]).all() for q in (p0, p) for key in fields))
    mass_exact = bool(np.array_equal(p0["mass"], p["mass"]))
    w0 = bool(np.all(p0["w2"] == 0) and np.all(p["w2"] == 0))
    positive = bool(all(np.all(q[key] > 0) for q in (p0, p)
                        for key in ("mass", "volume", "ie", "den", "pressure")))
    area0 = float(np.sum(p0["volume"].astype(float), dtype=np.longdouble))
    area1 = float(np.sum(p["volume"].astype(float), dtype=np.longdouble))
    checks = {
        "finite": finite, "fixed_mass_exact": mass_exact, "W0_exact": w0,
        "positive_state": positive,
        "periodic_area_tiles_box": bool(abs(area0 - 36.) < 1e-6 and
                                        abs(area1 - 36.) < 1e-5),
    }
    if not all(checks.values()):
        raise ValueError(f"Native provenance prerequisites failed: {checks}")
    x0 = p0["x"].astype(float) - 3.
    y0 = p0["y"].astype(float) - 3.
    x = p["x"].astype(float) - 3.
    y = p["y"].astype(float) - 3.
    r0, rt = np.hypot(x0, y0), np.hypot(x, y)
    mass = p0["mass"].astype(float)
    density_change = p["den"].astype(float) / p0["den"].astype(float) - 1.
    gm1 = float(np.float32(1.666667)) - 1.
    entropy_log = (np.log(p["ie"].astype(float) / p0["ie"].astype(float)) +
                   gm1 * np.log(p["volume"].astype(float) /
                                p0["volume"].astype(float)))
    vr = ((x * p["vx"].astype(float) + y * p["vy"].astype(float)) /
          np.maximum(rt, 1e-300))
    edge0 = np.minimum.reduce((p0["x"].astype(float), p0["y"].astype(float),
                               6. - p0["x"].astype(float),
                               6. - p0["y"].astype(float)))
    groups = [
        provenance("initial_outer_r_ge_2p75", r0 >= 2.75, r0, rt,
                   mass, density_change, entropy_log, vr),
        provenance("initial_seam_edge_lt_0p2", edge0 < .2, r0, rt,
                   mass, density_change, entropy_log, vr),
        provenance("initial_buffer_2p35_to_2p75", (r0 >= 2.35) & (r0 < 2.75),
                   r0, rt, mass, density_change, entropy_log, vr),
    ]
    disk = (r0 >= .5) & (r0 < 2.)
    gap = float(np.min(rt[r0 >= 2.75]) - np.max(rt[disk]))
    return {
        "initial": str(initial.resolve()), "initial_sha256": digest(initial),
        "endpoint": str(endpoint.resolve()), "endpoint_sha256": digest(endpoint),
        "time": t, "cells": len(p), "checks": checks,
        "initial_area_sum": area0, "endpoint_area_sum": area1,
        "initial_disk_cells": int(np.sum(disk)),
        "current_radial_gap_outer_origin_to_disk_origin": gap,
        "groups": groups,
        "scope": "Generator-centre provenance and fixed mass for W0 periodic Kepler. Crossings are not face-integrated fluid flux and absence of crossings cannot exclude transmitted pressure/torque.",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("initial", type=Path)
    ap.add_argument("endpoint", type=Path)
    ap.add_argument("--time", type=float, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = analyze(args.initial, args.endpoint, args.time)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"time": result["time"],
                      "radial_gap": result["current_radial_gap_outer_origin_to_disk_origin"],
                      "groups": result["groups"]}, indent=2))


if __name__ == "__main__":
    main()

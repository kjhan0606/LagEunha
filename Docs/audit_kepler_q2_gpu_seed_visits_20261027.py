"""Count linked-cell seed visits for one accepted native Kepler snapshot.

This measures bin occupancy, not actual half-plane clips or GPU instructions.
The raw-face CUDA tessellator scans the same seed bins once for directional
preclipping and once for the complete seed pass.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


def summarize(values):
    return {
        "sites": int(values.size),
        "total": int(np.sum(values, dtype=np.int64)),
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "p95": float(np.percentile(values, 95)),
        "p99": float(np.percentile(values, 99)),
        "maximum": int(np.max(values)),
        "over_512": int(np.sum(values > 512)),
    }


def counts_for_pad(prefix, ix, iy, pad, mx, my):
    xlo = np.maximum(ix - pad, 0)
    xhi = np.minimum(ix + pad + 1, mx)
    ylo = np.maximum(iy - pad, 0)
    yhi = np.minimum(iy + pad + 1, my)
    return (prefix[yhi, xhi] - prefix[ylo, xhi]
            - prefix[yhi, xlo] + prefix[ylo, xlo] - 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("snapshot", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    particles, _ = ordered(args.snapshot, 55351)
    if len(np.unique(particles["u4if"])) != len(particles):
        raise ValueError("Repeated particle IDs invalidate the self-only exclusion")
    mx = my = 98
    cellsize = .0625
    xmin = ymin = -.0625
    x = particles["x"].astype(float)
    y = particles["y"].astype(float)
    ix = np.clip(((x - xmin) / cellsize).astype(int), 0, mx - 1)
    iy = np.clip(((y - ymin) / cellsize).astype(int), 0, my - 1)
    occupancy = np.zeros((my, mx), dtype=np.int64)
    np.add.at(occupancy, (iy, ix), 1)
    prefix = np.pad(occupancy, ((1, 0), (1, 0))).cumsum(0).cumsum(1)
    r = np.hypot(x - 3., y - 3.)
    disk = (r >= .5) & (r < 2.)
    metrics = {}
    for pad in (1, 2):
        visits = counts_for_pad(prefix, ix, iy, pad, mx, my)
        metrics[f"pad_{pad}"] = {
            "all": summarize(visits),
            "disk": summarize(visits[disk]),
            "other": summarize(visits[~disk]),
        }
    result = {
        "snapshot": str(args.snapshot.resolve()),
        "snapshot_sha256": digest(args.snapshot),
        "snapshot_time": snapshot_time(args.snapshot, len(particles)),
        "grid": {"mx": mx, "my": my, "xmin": xmin, "ymin": ymin,
                 "cellsize": cellsize},
        "sites": len(particles),
        "unique_ids": True,
        "metrics": metrics,
        "scope": "Exact linked-cell occupancy and seed-bin visit counts; excludes adaptive bins, preclip early exits, clipping costs, and final wall planes."
    }
    with args.output.open("x") as out:
        json.dump(result, out, indent=2)
        out.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

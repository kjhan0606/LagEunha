"""Predict ring-coherent q2 W0 density change from exact circular site orbits."""

import hashlib
import json
from pathlib import Path
import struct
import time

import numpy as np

from audit_kepler_sparse_native_ballistic_20261027 import disk_geometry


DOCS = Path(__file__).resolve().parent
PREFIX = "KEPLER_Q2_BOUNDED_RINGCOHERENT005_R275_N1024_2026-10-27"
OUTPUT = DOCS / "KEPLER_Q2_BOUNDED_RINGCOHERENT005_BALLISTIC_T65_2026-10-27.json"
N = 55351


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def read_xy(raw):
    if len(raw) != 20 + 16*N or struct.unpack_from("<idd", raw) != (N, 6., 6.):
        raise ValueError("Unexpected ring glass header")
    x = np.frombuffer(raw, dtype="<f8", count=N, offset=20)
    y = np.frombuffer(raw, dtype="<f8", count=N, offset=20+8*N)
    return x, y


def read_area(raw):
    if len(raw) != 20 + 8*N or struct.unpack_from("<idd", raw) != (N, 6., 6.):
        raise ValueError("Unexpected ring area header")
    return np.frombuffer(raw, dtype="<f8", count=N, offset=20)


def metrics(area0, area, take):
    change = area0[take]/area[take]-1.
    weights = area0[take]
    return {
        "cells": int(np.sum(take)),
        "abs_relative_density_area_weighted": float(np.average(np.abs(change), weights=weights)),
        "abs_relative_density_p95": float(np.quantile(np.abs(change), .95)),
        "abs_relative_density_max": float(np.max(np.abs(change))),
    }


def main():
    started = time.monotonic()
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    manifest = json.loads((DOCS / f"{PREFIX}.json").read_text())
    glass = (DOCS / f"{PREFIX}.glass.bin").read_bytes()
    sidecar = (DOCS / f"{PREFIX}.areas.bin").read_bytes()
    if sha(glass) != manifest["glass_sha256"] or sha(sidecar) != manifest["areas_sha256"]:
        raise ValueError("Changed ring-coherent IC input")
    x, y = read_xy(glass)
    area0 = read_area(sidecar)
    r0 = np.hypot(x-3., y-3.)
    disk = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    if len(disk) != 49038 or not np.all(area0 > 0):
        raise ValueError("Unexpected ring-coherent disk")
    _, area, safety = disk_geometry({"x": x, "y": y}, 6.5, disk)
    if not np.all(area > 0):
        raise ValueError("Nonpositive exact-orbit disk area")
    result = {
        "glass_sha256": sha(glass), "areas_sha256": sha(sidecar),
        "time": 6.5, "disk": metrics(area0[disk], area, np.ones(len(disk), dtype=bool)),
        "inner": metrics(area0[disk], area, r0[disk] < .75),
        "minimum_periodic_image_safety_ratio": safety,
        "elapsed_seconds": time.monotonic()-started,
        "comparison_existing_ballistic_disk": {
            "aligned": .00018572657571010068,
            "independent_site_jitter005": .0014787874822204416,
        },
        "scope": "Exact softened circular site trajectories and independently rebuilt W0 disk Voronoi areas; no native force, entropy, or long-time stability conclusion.",
    }
    with OUTPUT.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

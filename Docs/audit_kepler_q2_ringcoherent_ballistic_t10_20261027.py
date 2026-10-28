"""Extend the ring-coherent exact-orbit geometry control to native t=10."""

import json
from pathlib import Path
import time

import numpy as np

from audit_kepler_q2_ringcoherent_ballistic_t65_20261027 import (
    DOCS, PREFIX, N, metrics, read_area, read_xy, sha,
)
from audit_kepler_sparse_native_ballistic_20261027 import disk_geometry


OUTPUT = DOCS / "KEPLER_Q2_BOUNDED_RINGCOHERENT005_BALLISTIC_T10_2026-10-27.json"


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
    _, area, safety = disk_geometry({"x": x, "y": y}, 10., disk)
    result = {
        "glass_sha256": sha(glass), "areas_sha256": sha(sidecar),
        "time": 10., "disk": metrics(area0[disk], area, np.ones(len(disk), dtype=bool)),
        "inner": metrics(area0[disk], area, r0[disk] < .75),
        "minimum_periodic_image_safety_ratio": safety,
        "elapsed_seconds": time.monotonic()-started,
        "comparison_existing_ballistic_disk": {
            "aligned": .0001942900523654185,
            "independent_site_jitter005": .001471921993474236,
        },
        "scope": "Exact softened circular site trajectories and independently rebuilt W0 disk Voronoi areas; no native force, entropy, or long-time stability conclusion.",
    }
    with OUTPUT.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

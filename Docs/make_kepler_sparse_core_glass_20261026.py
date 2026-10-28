"""Remove redundant low-mass central sites from the native q=2 glass IC.

Keep the disk and the exterior bitwise identical. Six sites on the outermost
ring below r=0.1 remain, evenly spaced in azimuth, so the center is covered
by gas-bearing Voronoi cells and no generator sits at r=0. This is an isolated
time-step control, not yet an equal-mass IC or a proven long-run solution.
"""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import cKDTree


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    raw = args.source.read_bytes()
    n, lx, ly = struct.unpack_from("<idd", raw)
    if (n, lx, ly) != (384*384, 6.0, 6.0) or len(raw) != 20+16*n:
        raise ValueError("Expected the complete native q=2 384² glass")
    x = np.frombuffer(raw, dtype="<f8", count=n, offset=20)
    y = np.frombuffer(raw, dtype="<f8", count=n, offset=20+8*n)
    r = np.hypot(x-3., y-3.)
    inner = r < .1
    first_ring = np.flatnonzero(inner & (abs(r-np.max(r[inner])) < 1e-9))
    if len(first_ring) < 6:
        raise ValueError("Cannot find six sites on outermost inner-core ring")
    angle = np.mod(np.arctan2(y[first_ring]-3., x[first_ring]-3.), 2*np.pi)
    wanted = np.arange(6)*2*np.pi/6
    chosen = []
    for target in wanted:
        j = int(np.argmin(abs(np.angle(np.exp(1j*(angle-target))))))
        chosen.append(int(first_ring[j]))
    if len(set(chosen)) != 6:
        raise ValueError("Repeated retained core site")
    keep = ~inner
    keep[chosen] = True
    xy = np.column_stack((x[keep], y[keep]))
    nearest = cKDTree(xy).query(xy, k=2)[0][:, 1]
    minimum_r = float(np.min(np.hypot(xy[:, 0]-3., xy[:, 1]-3.)))
    if not np.isfinite(nearest).all() or np.min(nearest) < .005 or minimum_r < .08:
        raise ValueError("Sparse glass has a collision or unexpectedly small radius")
    if np.sum((r>=.5) & (r<2.)) != np.sum(keep & (r>=.5) & (r<2.)):
        raise ValueError("Disk generators changed")
    with args.output.open("xb") as stream:
        stream.write(struct.pack("<idd", len(xy), lx, ly))
        stream.write(np.asarray(xy[:, 0], dtype="<f8").tobytes())
        stream.write(np.asarray(xy[:, 1], dtype="<f8").tobytes())
    result = {
        "source": str(args.source),
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "output": str(args.output),
        "output_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
        "source_cells": n,
        "output_cells": len(xy),
        "removed_inner_cells": int(n-len(xy)),
        "retained_inner_cells": 6,
        "retained_core_site_indices": chosen,
        "minimum_radius": minimum_r,
        "minimum_pair_distance": float(np.min(nearest)),
        "disk_generators_unchanged": True,
        "requires_native_env": "EUNHA_ALLOW_SPARSE_GLASS=1 and EUNHA_KEPLER_MASS_FROM_CELL_VOLUME=1",
        "scope": "Core-time-step control; equal-mass and native stability not established",
    }
    with args.output.with_suffix(".json").open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

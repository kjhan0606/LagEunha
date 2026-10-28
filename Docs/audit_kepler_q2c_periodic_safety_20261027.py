"""Certify q2c disk polygons against omitted periodic images at early times."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from audit_kepler_polar_long_geometry_20261026 import initial, one_time


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    p, geometry = initial(args.snapshot)
    if len(p) != 384 * 384:
        raise ValueError("Expected complete q2c 384x384 snapshot")
    if not all(np.isfinite(p[key]).all() for key in ("mass", "volume", "den")):
        raise ValueError("Nonfinite initial mass, volume or density")
    if not all(np.all(p[key] > 0) for key in ("mass", "volume", "den")):
        raise ValueError("Nonpositive initial mass, volume or density")
    area_sum = float(np.sum(p["volume"].astype(float), dtype=np.longdouble))
    rho_error = float(np.max(abs(p["den"].astype(float) /
                                 (p["mass"].astype(float) /
                                  p["volume"].astype(float)) - 1.)))
    if abs(area_sum - 36.) > 1e-5 or rho_error >= 1e-6:
        raise ValueError("Initial periodic area or density relation failed")
    result = {
        "initial": str(args.snapshot),
        "initial_sha256": sha256(args.snapshot),
        "initial_area_sum": area_sum,
        "rho_mass_over_volume_max_relative": rho_error,
        "times": [one_time(p, geometry, t) for t in (0., .825)],
        "scope": "Exact-orbit W0 q2c disk geometry with omitted-periodic-image safety certificate at t=0 and .825; no native long-time fluid stability claim",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"minimum_safety_ratios": [
        row["minimum_periodic_image_safety_ratio"] for row in result["times"]
    ]}))


if __name__ == "__main__":
    main()

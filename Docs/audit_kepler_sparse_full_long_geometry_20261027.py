"""Same-time exact-orbit Voronoi geometry of full and corrected sparse disks.

This uses certified native t=0 areas and fixed masses. It moves generators on
the analytic softened circular orbits but does not evolve hydrodynamics.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from audit_kepler_polar_long_geometry_20261026 import one_time
from plot_kepler_native_density_long_20261026 import complete_ordered


TIMES = (0., .825, 2.221441469079183, 6.5, 10., 120., 600.)
COUNTS = {"full_q2": 384 * 384, "sparse_narrow": 54312, "sparse_wide": 55931}


def digest(path):
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def evaluate(name, path):
    p, meta = complete_ordered(path)
    if len(p) != COUNTS[name] or float(meta["t"]) != 0.:
        raise ValueError(f"{name}: expected exact complete native t=0 snapshot")
    if not np.all(p["w2"] == 0):
        raise ValueError(f"{name}: expected W=0 Voronoi geometry")
    if not all(np.isfinite(p[key]).all() for key in ("x", "y", "mass", "volume", "den")):
        raise ValueError(f"{name}: nonfinite initial state")
    if not (np.all(p["mass"] > 0) and np.all(p["volume"] > 0) and np.all(p["den"] > 0)):
        raise ValueError(f"{name}: nonpositive initial state")
    area_sum = float(np.sum(p["volume"].astype(float), dtype=np.longdouble))
    if abs(area_sum - 36.) > 1e-5:
        raise ValueError(f"{name}: native t=0 areas do not tile 6x6 box")
    rho_error = float(np.max(abs(p["den"].astype(float) /
                                 (p["mass"].astype(float) / p["volume"].astype(float)) - 1.)))
    if rho_error > 1e-6:
        raise ValueError(f"{name}: native t=0 density differs from mass/volume")
    x, y = p["x"].astype(float) - 3., p["y"].astype(float) - 3.
    radius = np.hypot(x, y)
    selected = np.flatnonzero((radius >= .5) & (radius < 2.))
    geometry = (x, y, radius, selected)
    return {"initial": str(path), "sha256": digest(path), "sites": len(p),
            "initial_area_sum": area_sum, "rho_mass_volume_max_relative": rho_error,
            "disk_sites": len(selected),
            "times": [one_time(p, geometry, t) for t in TIMES]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("full_q2", type=Path)
    parser.add_argument("sparse_narrow", type=Path)
    parser.add_argument("sparse_wide", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    cases = {name: evaluate(name, getattr(args, name)) for name in COUNTS}
    comparison = []
    for k, t in enumerate(TIMES):
        full = cases["full_q2"]["times"][k]["bins"][-1]["density_change_abs_mass_weighted"]
        row = {"time": t, "full_q2_disk_density_abs_mass_weighted": full}
        for name in ("sparse_narrow", "sparse_wide"):
            sparse = cases[name]["times"][k]["bins"][-1]["density_change_abs_mass_weighted"]
            row[name + "_disk_density_abs_mass_weighted"] = sparse
            row[name + "_over_full_q2"] = sparse / full if full else None
        comparison.append(row)
    result = {
        "cases": cases, "comparison": comparison,
        "scope": "Exact softened circular-orbit Voronoi areas of native W=0 fixed-mass initial states; disk polygons certified independent of omitted periodic images. No native gas, pressure, entropy, boundary evolution, or t=600 stability conclusion.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"comparison": comparison}, indent=2))


if __name__ == "__main__":
    main()

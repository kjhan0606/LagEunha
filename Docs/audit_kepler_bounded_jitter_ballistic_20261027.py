"""Exact-orbit Voronoi geometry for aligned and 5% jitter bounded q2 ICs."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_native_ballistic_20261027 import disk_geometry
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N = 55351
TIMES = (0.825, 6.5)


def calculate(path):
    if abs(snapshot_time(path, N)) > 1e-12:
        raise ValueError(f"Expected complete t=0 native state: {path}")
    state, _ = ordered(path, N)
    r = np.hypot(state["x"].astype(float)-3, state["y"].astype(float)-3)
    selected = np.flatnonzero((r >= .5) & (r < 2.))
    if len(selected) != 49038 or np.any(state["w2"] != 0):
        raise ValueError("Unexpected disk population or nonzero weight")
    initial_area = state["volume"][selected].astype(float)
    mass = state["mass"][selected].astype(float)
    metrics = {}
    for time in TIMES:
        _, area, safety = disk_geometry(state, time, selected)
        relative = initial_area / area - 1
        metrics[str(time)] = {
            "disk_abs_relative_density_mass_weighted": float(
                np.sum(np.abs(relative) * mass, dtype=np.longdouble)
                / np.sum(mass, dtype=np.longdouble)),
            "disk_p95_abs_relative_density": float(np.quantile(np.abs(relative), .95)),
            "disk_max_abs_relative_density": float(np.max(np.abs(relative))),
            "minimum_periodic_image_safety_ratio": float(safety),
        }
    return {"snapshot": str(path.resolve()), "sha256": digest(path),
            "disk_sites": len(selected), "metrics": metrics}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("aligned", type=Path)
    ap.add_argument("jitter", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    aligned = calculate(args.aligned)
    jitter = calculate(args.jitter)
    ratio = {}
    for time in TIMES:
        key = str(time)
        ratio[key] = (jitter["metrics"][key]["disk_abs_relative_density_mass_weighted"]
                      / aligned["metrics"][key]["disk_abs_relative_density_mass_weighted"])
    result = {
        "aligned": aligned, "jitter": jitter,
        "jitter_to_aligned_disk_density_abs_ratio": ratio,
        "scope": "Ballistic generator trajectories and independently rebuilt Voronoi geometry only; no native force, entropy, or long-time stability conclusion",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

"""Trace selected q2 low-CFL cell IDs from t=.001 to t=2.0.

This is a nearest-real-generator diagnostic, not a reconstruction of the full
periodic Voronoi face graph. The saved per-cell dt may have come from an RK
stage geometry, so the nearest-pair dt estimate is only a consistency check.
"""

import argparse
import json
import struct
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


N_EXPECTED = 384 * 384
STRIDE_EXPECTED = 576
TARGET_IDS = (8512, 146887, 10127, 143520, 76397, 127494)
GAMMA = 1.666667
COURANT = 0.4
ETA = 0.06
EPS = 0.01175
CX = 3.0
CY = 3.0

DTYPE = np.dtype({
    "names": ["u4if", "x", "y", "mass", "vx", "vy", "den", "pressure", "ie", "volume", "dt"],
    "formats": ["<u8", "<f8", "<f8", "<f8", "<f8", "<f8", "<f4", "<f4", "<f4", "<f4", "<f4"],
    "offsets": [8, 24, 32, 48, 56, 64, 104, 108, 116, 132, 128],
    "itemsize": STRIDE_EXPECTED,
})


def load(path: Path):
    with path.open("rb") as stream:
        header = stream.read(16)
        if len(header) != 16:
            raise ValueError(f"Truncated native header: {path}")
        magic, n, stride, stress = struct.unpack("<4i", header)
        if (magic, n, stride, stress) != (-1, N_EXPECTED, STRIDE_EXPECTED, 1):
            raise ValueError(f"Unexpected native layout: {(magic, n, stride, stress)}")
        stream.seek(16 + n * stride)
        footer = stream.read(52)
        if len(footer) != 52:
            raise ValueError(f"Missing native footer: {path}")
        time = struct.unpack_from("<d", footer)[0]
    if path.stat().st_size != 16 + n * stride + 52:
        raise ValueError(f"Incomplete native snapshot: {path}")

    data = np.memmap(path, mode="r", dtype=DTYPE, offset=16, shape=(n,))
    ids = np.asarray(data["u4if"], dtype=np.uint64) & np.uint64((1 << 56) - 1)
    if not np.array_equal(np.sort(ids), np.arange(n, dtype=np.uint64)):
        raise ValueError(f"Cell IDs are not unique and complete: {path}")
    if not np.isfinite(np.asarray(data["dt"], dtype=float)).all():
        raise ValueError(f"Nonfinite saved hydro timestep: {path}")
    return data, ids.astype(np.int64), time


def inspect(path: Path) -> dict:
    data, ids, time = load(path)
    row_by_id = np.empty(N_EXPECTED, dtype=np.int64)
    row_by_id[ids] = np.arange(N_EXPECTED)
    xy = np.column_stack((data["x"], data["y"]))
    velocity = np.column_stack((data["vx"], data["vy"]))
    sound_speed = np.sqrt(
        GAMMA * np.asarray(data["pressure"], dtype=float)
        / np.asarray(data["den"], dtype=float)
    )
    tree = cKDTree(xy)
    rows = row_by_id[np.asarray(TARGET_IDS)]
    distances, neighbors = tree.query(xy[rows], k=2)
    observations = []

    for index, row in enumerate(rows):
        ident = TARGET_IDS[index]
        neighbor = int(neighbors[index, 1])
        delta = xy[neighbor] - xy[row]
        gap = float(distances[index, 1])
        normal = delta / gap
        radial = xy[row] - np.array([CX, CY])
        radius = float(np.linalg.norm(radial))
        radial_unit = radial / radius
        azimuthal_unit = np.array([-radial_unit[1], radial_unit[0]])
        circular_speed = radius / (radius * radius + EPS * EPS) ** 0.75
        vr = float(velocity[row] @ radial_unit)
        vphi = float(velocity[row] @ azimuthal_unit)

        neighbor_radius = float(np.linalg.norm(xy[neighbor] - np.array([CX, CY])))
        relative_normal_velocity = float((velocity[neighbor] - velocity[row]) @ normal)
        signal_speed = float(
            sound_speed[row] + sound_speed[neighbor]
            - min(0.0, relative_normal_velocity)
        )
        nearest_pair_dt = 2.0 * COURANT * gap / signal_speed

        observations.append({
            "id": ident,
            "radius": radius,
            "vr_over_vc": vr / circular_speed,
            "vphi_over_vc": vphi / circular_speed,
            "saved_hydro_dt": float(data["dt"][row]),
            "volume": float(data["volume"][row]),
            "mass": float(data["mass"][row]),
            "density": float(data["den"][row]),
            "pressure": float(data["pressure"][row]),
            "internal_energy": float(data["ie"][row]),
            "nearest_real_neighbor_id": int(ids[neighbor]),
            "nearest_gap": gap,
            "nearest_neighbor_radius": neighbor_radius,
            "radial_gap": neighbor_radius - radius,
            "relative_normal_velocity": relative_normal_velocity,
            "nearest_pair_dt_estimate": nearest_pair_dt,
        })

    return {"snapshot": str(path), "time": time, "cells": N_EXPECTED, "observations": observations}


def compare(initial: dict, evolved: dict) -> list[dict]:
    old = {item["id"]: item for item in initial["observations"]}
    new = {item["id"]: item for item in evolved["observations"]}
    results = []
    for ident in TARGET_IDS:
        a, b = old[ident], new[ident]
        entropy_a = a["pressure"] / a["density"] ** GAMMA
        entropy_b = b["pressure"] / b["density"] ** GAMMA
        results.append({
            "id": ident,
            "initial": a,
            "evolved": b,
            "ratios": {
                "mass": b["mass"] / a["mass"],
                "volume": b["volume"] / a["volume"],
                "density": b["density"] / a["density"],
                "pressure": b["pressure"] / a["pressure"],
                "internal_energy": b["internal_energy"] / a["internal_energy"],
                "entropy_proxy": entropy_b / entropy_a,
            },
        })
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("evolved", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    initial = inspect(args.initial)
    evolved_snapshots = [inspect(path) for path in args.evolved]
    if any(not item["time"] > initial["time"] for item in evolved_snapshots):
        raise ValueError("Every evolved snapshot must be later than the initial snapshot")
    if any(a["time"] >= b["time"] for a, b in zip(evolved_snapshots, evolved_snapshots[1:])):
        raise ValueError("Evolved snapshots must be supplied in strictly increasing time order")
    result = {
        "scope": (
            "Selected-cell lineage and nearest real-generator diagnostic. It does not "
            "reconstruct periodic images or the exact RK-stage Voronoi face graph; "
            "nearest-pair dt estimates are not exact provenance for saved per-cell dt."
        ),
        "formula": "dt_pair=2*C*distance/(cs_i+cs_j-min(0,dv_dot_e)); cs=sqrt(gamma*P/rho)",
        "parameters": {"gamma": GAMMA, "Courant": COURANT, "eta": ETA, "eps": EPS, "center": [CX, CY]},
        "initial_snapshot": initial,
        "evolved_snapshots": evolved_snapshots,
        "lineage": [
            {"time": item["time"], "comparisons": compare(initial, item)}
            for item in evolved_snapshots
        ],
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

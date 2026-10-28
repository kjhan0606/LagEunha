"""Check native CUDA pressure work and geometry on fixed Kepler wall faces."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N = 55351
SIDES = 1024
RADIUS = 2.75
APOTHEM = RADIUS * np.cos(np.pi / SIDES)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    state_path = args.run / "khout.000000.dat"
    faces_path = args.run / "face_dump_r0_call1.tsv"
    if abs(snapshot_time(state_path, N)) > 1e-12:
        raise ValueError("Missing complete native initial state")
    state, _ = ordered(state_path, N)
    faces = np.genfromtxt(faces_path, names=True, delimiter="\t", dtype=None,
                          encoding="utf-8")
    faces = np.atleast_1d(faces)
    if len(faces) == 0:
        raise ValueError("No physical wall faces")

    ids = faces["id_i"].astype(int)
    sides = -1000 - faces["j"].astype(int)
    nx = np.cos((sides + 0.5) * 2 * np.pi / SIDES)
    ny = np.sin((sides + 0.5) * 2 * np.pi / SIDES)
    dsx = faces["dSx"].astype(float)
    dsy = faces["dSy"].astype(float)
    length = np.hypot(dsx, dsy)
    xmid = 0.5 * (faces["x1"].astype(float) + faces["x2"].astype(float)) - 3
    ymid = 0.5 * (faces["y1"].astype(float) + faces["y2"].astype(float)) - 3
    vx = state["vx"][ids].astype(float)
    vy = state["vy"][ids].astype(float)
    pressure = state["pressure"][ids].astype(float)
    pstar = faces["pstar"].astype(float)
    sweep = faces["udotds"].astype(float)
    work = faces["work"].astype(float)

    endpoint_plane_error = np.maximum(
        np.abs((faces["x1"].astype(float) - 3) * nx
               + (faces["y1"].astype(float) - 3) * ny - APOTHEM),
        np.abs((faces["x2"].astype(float) - 3) * nx
               + (faces["y2"].astype(float) - 3) * ny - APOTHEM))
    outward_alignment_error = np.hypot(dsx / length - nx, dsy / length - ny)
    sweep_error = sweep + vx * dsx + vy * dsy
    work_error = work + pstar * sweep
    p_rel_error = np.abs(pstar / pressure - 1)
    wall_force = np.array([-np.sum(pstar * dsx, dtype=np.longdouble),
                           -np.sum(pstar * dsy, dtype=np.longdouble)], dtype=float)
    wall_torque = -float(np.sum(pstar * (xmid * dsy - ymid * dsx),
                                dtype=np.longdouble))
    total_pressure_scale = float(np.sum(np.abs(pstar) * length,
                                        dtype=np.longdouble))
    torque_scale = float(np.sum(np.abs(pstar) * length * np.hypot(xmid, ymid),
                                  dtype=np.longdouble))

    checks = {
        "wall_rows_only": bool(np.all(faces["j"] <= -1000)
                               and np.all(faces["wall"] == 2)),
        "valid_side_ids": bool(np.all((sides >= 0) & (sides < SIDES))),
        "positive_face_lengths": bool(np.all(np.isfinite(length) & (length > 0))),
        "faces_on_polygon": bool(np.max(endpoint_plane_error) < 1e-10),
        "outward_area_vectors": bool(np.max(outward_alignment_error) < 1e-10),
        "native_relative_sweep": bool(np.max(np.abs(sweep_error)) < 1e-12),
        "native_pressure_work": bool(np.max(np.abs(work_error)) < 1e-16),
        "native_cell_pressure": bool(np.max(p_rel_error) < 1e-6),
        "uniform_wall_net_force": bool(np.linalg.norm(wall_force)
                                       / total_pressure_scale < 1e-6),
        "uniform_wall_net_torque": bool(abs(wall_torque)
                                        / torque_scale < 1e-6),
    }
    result = {
        "run": str(args.run.resolve()),
        "snapshot_sha256": digest(state_path),
        "faces_sha256": digest(faces_path),
        "wall_face_rows": len(faces),
        "wall_cells": len(np.unique(ids)),
        "polygon_sides_covered": len(np.unique(sides)),
        "max_endpoint_plane_error": float(np.max(endpoint_plane_error)),
        "max_outward_alignment_error": float(np.max(outward_alignment_error)),
        "max_relative_sweep_error": float(np.max(np.abs(sweep_error))),
        "max_pressure_work_error": float(np.max(np.abs(work_error))),
        "max_pressure_relative_error": float(np.max(p_rel_error)),
        "net_wall_pressure_force": wall_force.tolist(),
        "net_wall_pressure_torque": wall_torque,
        "total_wall_work_rate": float(np.sum(work, dtype=np.longdouble)),
        "checks": checks,
        "scope": "Initial native CUDA wall faces and work only; no generator-wall collision",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))
    if not all(checks.values()):
        raise SystemExit("Native wall-face audit failed")


if __name__ == "__main__":
    main()

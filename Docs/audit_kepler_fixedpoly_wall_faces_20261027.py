"""Audit one native CUDA fixed-polygon wall-face dump.

The input is a compact GFS_GPU_FACE_DUMP_WALL_ONLY=1 capture made with
KEPLER_FIXEDPOLY_WALL_DUMP_2026-10-27.patch.  Its `work` column is the
relative-sweep correction in the cell energy equation, not work done by a
stationary wall on the fluid.
"""

import argparse
import csv
import json
import math
from pathlib import Path


SIDES = 1024
RADIUS = 2.75
CENTER = (3.0, 3.0)
APOTHEM = RADIUS * math.cos(math.pi / SIDES)
EDGE_LENGTH = 2.0 * RADIUS * math.sin(math.pi / SIDES)
PERIMETER = SIDES * EDGE_LENGTH


def audit(path: Path, initial_uniform: bool = True):
    side_lengths = [0.0] * SIDES
    side_counts = [0] * SIDES
    length_total = 0.0
    max_plane_error = 0.0
    max_normal_error = 0.0
    max_work_error = 0.0
    max_pressure_relative = 0.0
    fx = fy = torque = work_signed = work_absolute = 0.0
    sweep_signed = sweep_absolute = 0.0
    pressure_reference = None
    calls = set()
    owners = set()
    rows = 0
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"call", "j", "id_i", "id_j", "wall", "pstar",
                    "udotds", "work", "dSx", "dSy", "x1", "y1", "x2", "y2"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError("Missing wall-face dump columns")
        for row in reader:
            rows += 1
            calls.add(int(row["call"]))
            owner = int(row["id_i"])
            owners.add(owner)
            j = int(row["j"])
            if int(row["wall"]) != 2 or int(row["id_j"]) != j:
                raise ValueError(f"Nonphysical or mismatched wall face at row {rows}")
            side = -1000 - j
            if not 0 <= side < SIDES:
                raise ValueError(f"Invalid polygon support {j} at row {rows}")
            pressure = float(row["pstar"])
            sweep = float(row["udotds"])
            work = float(row["work"])
            dsx, dsy = float(row["dSx"]), float(row["dSy"])
            x1, y1 = float(row["x1"]), float(row["y1"])
            x2, y2 = float(row["x2"]), float(row["y2"])
            values = (pressure, sweep, work, dsx, dsy, x1, y1, x2, y2)
            if not all(math.isfinite(value) for value in values) or pressure <= 0:
                raise ValueError(f"Nonfinite or nonpositive wall state at row {rows}")
            edge_length = math.hypot(dsx, dsy)
            if edge_length <= 0:
                raise ValueError(f"Zero wall face at row {rows}")
            if abs(dsx - (y2 - y1)) > 1e-11 or abs(dsy + (x2 - x1)) > 1e-11:
                raise ValueError(f"Inconsistent face vector at row {rows}")
            angle = (side + 0.5) * 2.0 * math.pi / SIDES
            nx, ny = math.cos(angle), math.sin(angle)
            mx, my = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
            plane_error = abs(nx * (mx - CENTER[0]) + ny * (my - CENTER[1]) - APOTHEM)
            normal_error = math.hypot(dsx / edge_length - nx, dsy / edge_length - ny)
            max_plane_error = max(max_plane_error, plane_error)
            max_normal_error = max(max_normal_error, normal_error)
            max_work_error = max(max_work_error, abs(work + pressure * sweep))
            if pressure_reference is None:
                pressure_reference = pressure
            max_pressure_relative = max(max_pressure_relative,
                                        abs(pressure / pressure_reference - 1.0))
            side_lengths[side] += edge_length
            side_counts[side] += 1
            length_total += edge_length
            fx -= pressure * dsx
            fy -= pressure * dsy
            torque -= pressure * ((mx - CENTER[0]) * dsy - (my - CENTER[1]) * dsx)
            work_signed += work
            work_absolute += abs(work)
            sweep_signed += sweep
            sweep_absolute += abs(sweep)
    if rows == 0:
        raise ValueError("No fixed-polygon wall faces in dump")
    max_side_relative_error = max(abs(length / EDGE_LENGTH - 1.0)
                                  for length in side_lengths)
    traction_scale = pressure_reference * PERIMETER
    torque_scale = traction_scale * RADIUS
    checks = {
        "single_force_call": len(calls) == 1,
        "all_polygon_sides_present": all(side_counts),
        "wall_perimeter_partition": abs(length_total / PERIMETER - 1.0) < 1e-7,
        "each_wall_side_partition": max_side_relative_error < 1e-6,
        "wall_faces_lie_on_support": max_plane_error < 1e-8,
        "wall_face_normals_outward": max_normal_error < 1e-8,
        "relative_sweep_energy_term": max_work_error < 1e-11 * traction_scale,
    }
    if initial_uniform:
        checks.update({
            "initial_pressure_uniform": max_pressure_relative < 1e-5,
            "initial_wall_force_cancels": math.hypot(fx, fy) < 1e-5 * traction_scale,
            "initial_wall_torque_cancels": abs(torque) < 1e-5 * torque_scale,
        })
    return {
        "input": str(path.resolve()), "force_calls": sorted(calls),
        "wall_faces": rows, "wall_cells": len(owners),
        "polygon_sides": SIDES, "sides_touched": sum(bool(n) for n in side_counts),
        "expected_perimeter": PERIMETER, "observed_perimeter": length_total,
        "maximum_side_length_relative_error": max_side_relative_error,
        "maximum_support_plane_error": max_plane_error,
        "maximum_outward_normal_error": max_normal_error,
        "pressure_reference": pressure_reference,
        "maximum_pressure_relative_change": max_pressure_relative,
        "wall_force_on_fluid": [fx, fy], "wall_torque_on_fluid": torque,
        "relative_sweep_sum": sweep_signed,
        "relative_sweep_absolute_sum": sweep_absolute,
        "energy_correction_sum": work_signed,
        "energy_correction_absolute_sum": work_absolute,
        "maximum_energy_correction_formula_error": max_work_error,
        "checks": checks,
        "scope": "One CUDA force call; static-wall traction and geometric energy correction. Physical wall power is zero.",
        "requires_uniform_initial_pressure": initial_uniform,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("face_dump", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-nonuniform", action="store_true",
                        help="check geometry and work, report nonzero wall force/torque at later times")
    args = parser.parse_args()
    receipt = audit(args.face_dump, initial_uniform=not args.allow_nonuniform)
    with args.output.open("x") as stream:
        json.dump(receipt, stream, indent=2)
        stream.write("\n")
    print(json.dumps(receipt, indent=2))
    if not all(receipt["checks"].values()):
        raise SystemExit("Fixed-polygon wall-face audit failed")


if __name__ == "__main__":
    main()

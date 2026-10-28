#!/usr/bin/env python3
"""Summarize intrinsic spin and kinetic-energy partitions in a FullMoment state.

This is a diagnostic only. It reads accepted checkpoint arrays and does not
modify them or alter the physical model.  The state layout is the 20-double
per-cell layout in Exam/full_moment_state_io.c and
Exam/full_moment_coupled_state.c.
"""

from __future__ import annotations

import argparse
import array
import json
import math
import struct
import sys
from pathlib import Path

MAGIC = 0x31544154534D464C


def read_array(path: Path, typecode: str, count: int) -> array.array:
    values = array.array(typecode)
    with path.open("rb") as stream:
        values.fromfile(stream, count)
        if stream.read(1):
            raise ValueError(f"unexpected trailing bytes in {path}")
    if values.itemsize != 8:
        raise RuntimeError(f"expected 8-byte array elements, got {values.itemsize}: {path}")
    if sys.byteorder != "little":
        values.byteswap()
    return values


def summarize(checkpoint: Path) -> dict[str, object]:
    checkpoint = checkpoint.resolve()
    raw_header = (checkpoint / "header.bin").read_bytes()
    if len(raw_header) != 128:
        raise ValueError(f"invalid header length in {checkpoint}: {len(raw_header)}")
    header = struct.unpack("<8Qd7Q", raw_header)
    magic, version, cells, epoch, token, state_elements, ld_bytes, ld_mantissa = header[:8]
    time = header[8]
    if magic != MAGIC or version != 1 or cells < 1 or state_elements != 20 * cells + 4:
        raise ValueError(f"invalid FullMoment state header: {checkpoint}")

    mass = read_array(checkpoint / "mass.bin", "d", cells)
    state = read_array(checkpoint / "state.bin", "d", state_elements)
    ids = read_array(checkpoint / "ids.bin", "q", cells)

    names = (
        "mass", "px", "py", "angular_momentum_z", "total_energy",
        "kinetic_second_moment", "bulk_kinetic", "internal_kinetic",
        "affine_kinetic", "rigid_rotation_energy", "non_affine_kinetic",
    )
    totals = {name: [] for name in names}
    min_inertia = math.inf
    min_thermal = math.inf
    max_abs_omega = -1.0
    max_abs_omega_cell = None

    for i, m in enumerate(mass):
        a = state[20 * i : 20 * i + 20]
        if not math.isfinite(m) or m <= 0.0 or not all(math.isfinite(x) for x in a):
            raise ValueError(f"non-finite or nonpositive state at cell index {i}")
        if a[6] != a[7] or a[16] != a[17]:
            raise ValueError(f"asymmetric raw moments at cell index {i}")

        qx, qy = a[3], a[4]
        px, py = a[9], a[10]
        sxx = a[5] - qx * qx / m
        sxy = a[6] - qx * qy / m
        syy = a[8] - qy * qy / m
        cxx = a[11] - qx * px / m
        cxy = a[12] - qx * py / m
        cyx = a[13] - qy * px / m
        cyy = a[14] - qy * py / m
        inertia = sxx + syy
        determinant = sxx * syy - sxy * sxy
        if sxx <= 0.0 or determinant <= 0.0 or inertia <= 0.0:
            raise ValueError(f"nonpositive centered position covariance at cell index {i}")

        # L_int = C_xy - C_yx; K_rot is its rigid-rotation projection.
        l_internal = cxy - cyx
        omega = l_internal / inertia
        k_rotation = 0.5 * l_internal * l_internal / inertia

        # K_affine = 1/2 tr(C^T S^{-1} C), with C_ij = integral x_i v_j.
        inv00 = syy / determinant
        inv01 = -sxy / determinant
        inv11 = sxx / determinant
        affine_trace = 0.0
        for c0, c1 in ((cxx, cyx), (cxy, cyy)):
            affine_trace += c0 * (inv00 * c0 + inv01 * c1)
            affine_trace += c1 * (inv01 * c0 + inv11 * c1)
        k_affine = 0.5 * affine_trace

        k_internal = 0.5 * ((a[15] - px * px / m) + (a[18] - py * py / m))
        k_non_affine = k_internal - k_affine
        thermal = a[19] - 0.5 * (a[15] + a[18])
        row = (
            m, px, py, a[12] - a[13], a[19],
            0.5 * (a[15] + a[18]),
            0.5 * (px * px + py * py) / m,
            k_internal, k_affine, k_rotation, k_non_affine,
        )
        if not all(math.isfinite(x) for x in row) or not math.isfinite(thermal):
            raise ValueError(f"non-finite diagnostic at cell index {i}")
        for name, value in zip(names, row):
            totals[name].append(value)
        min_inertia = min(min_inertia, inertia)
        min_thermal = min(min_thermal, thermal)
        if abs(omega) > max_abs_omega:
            max_abs_omega = abs(omega)
            max_abs_omega_cell = int(ids[i])

    sums = {name: math.fsum(values) for name, values in totals.items()}
    internal = sums["internal_kinetic"]
    summary: dict[str, object] = {
        "checkpoint": str(checkpoint),
        "cells": cells,
        "time": time,
        "epoch": epoch,
        "run_token": token,
        "ld_bytes": ld_bytes,
        "ld_mantissa": ld_mantissa,
        **sums,
        "rigid_rotation_fraction_of_internal_kinetic": (
            sums["rigid_rotation_energy"] / internal if internal != 0.0 else None
        ),
        "max_abs_intrinsic_spin_omega": max_abs_omega,
        "max_abs_intrinsic_spin_omega_cell_id": max_abs_omega_cell,
        "minimum_cell_inertia": min_inertia,
        "minimum_cell_thermal": min_thermal,
        "kinetic_partition_residual": (
            internal - sums["affine_kinetic"] - sums["non_affine_kinetic"]
        ),
        "interpretation": (
            "Rotation is diagnosed from the existing off-diagonal raw moments. "
            "K_rot is a projection of existing kinetic energy, not an added energy."
        ),
    }
    accepted = checkpoint.parent / f"accepted_{checkpoint.name.rsplit('_', 1)[-1]}.json"
    if accepted.is_file():
        record = json.loads(accepted.read_text())
        summary["accepted_step"] = record.get("step")
        summary["step_wall_seconds"] = record.get("step_wall_seconds")
        summary["long_time_stability_verified"] = record.get("long_time_stability_verified")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoints", nargs="+", type=Path)
    args = parser.parse_args()
    result = [summarize(path) for path in args.checkpoints]
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

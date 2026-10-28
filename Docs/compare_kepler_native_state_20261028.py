#!/usr/bin/env python3
"""Compare two native RK4 Kepler snapshots while ignoring serialized pointers.

The native particle record includes linked-list pointers and ABI padding, so
raw file comparison is not a physical-state comparison. This reader models the
current XYZDBL + stress record layout and compares every stored numeric field.
"""

import argparse
import ctypes as C
import json
from pathlib import Path
import struct

import numpy as np


class NativeParticle(C.Structure):
    pass


_fields = [
    ("next", C.c_void_p),
    ("u4if", C.c_uint64),
    ("indt", C.c_char * 4),
]
_fields += [
    (name, C.c_double)
    for name in ("x", "y", "z", "mass", "vx", "vy", "vz", "ax", "ay", "az")
]
_fields += [
    (name, C.c_float)
    for name in (
        "den", "pressure", "avgNeighboringPressure", "ie", "die", "csound",
        "dt", "volume", "w2", "w2hydro", "w2old", "w2ceil",
    )
]
_rk4 = (
    "k1x", "k1y", "k1vx", "k1vy", "k2x", "k2y", "k2vx", "k2vy",
    "k3x", "k3y", "k3vx", "k3vy", "k4x", "k4y", "k4vx", "k4vy",
    "k1ie", "k2ie", "k3ie", "k4ie", "k1K", "k2K", "k3K", "k4K",
)
_fields += [(name, C.c_double) for name in _rk4]
_fields += [("w2backup", C.c_float)]
_stress = (
    "gUxx", "gUxy", "gUyx", "gUyy", "divv", "dPdx", "dPdy", "dRhodx",
    "dRhody", "tauxx", "tauxy", "tauyy", "divv_old", "alpha_cd",
    "vsig_max", "E_inv_xx", "E_inv_xy", "E_inv_yx", "E_inv_yy", "h_mfm",
    "mass_flux_accum", "v_smooth_x", "v_smooth_y", "lap_vx", "lap_vy", "K", "dK",
)
_fields += [(f"stress_{name}", C.c_double) for name in _stress]
_fields += [("bp", C.c_void_p)]
NativeParticle._fields_ = _fields

RECORD_SIZE = C.sizeof(NativeParticle)
PRIMARY_STATE = (
    "x", "y", "z", "mass", "vx", "vy", "vz", "den", "pressure",
    "avgNeighboringPressure", "ie", "die", "csound", "dt", "volume",
    "w2", "w2hydro", "w2old", "w2ceil",
)


def record_dtype():
    names, formats, offsets = [], [], []
    for name, field_type in NativeParticle._fields_:
        names.append(name)
        offsets.append(getattr(NativeParticle, name).offset)
        if name in ("next", "bp"):
            formats.append("<u8")
        elif name == "u4if":
            formats.append("<u8")
        elif name == "indt":
            formats.append("V4")
        elif field_type is C.c_double:
            formats.append("<f8")
        elif field_type is C.c_float:
            formats.append("<f4")
        else:
            raise TypeError(f"Unsupported native field {name}: {field_type}")
    return np.dtype({
        "names": names,
        "formats": formats,
        "offsets": offsets,
        "itemsize": RECORD_SIZE,
    })


def load_snapshot(path):
    path = Path(path)
    with path.open("rb") as stream:
        header = struct.unpack("<4i", stream.read(16))
    magic, count, stride, has_stress = header
    if (magic, stride, has_stress) != (-1, RECORD_SIZE, 1):
        raise ValueError(f"Unsupported native header {header} in {path}")
    if path.stat().st_size != 16 + count * stride + 52:
        raise ValueError(f"Incomplete native snapshot {path}")
    with path.open("rb") as stream:
        stream.seek(16 + count * stride)
        raw_meta = stream.read(52)
    t, dt, dtold, nx, ny, lx, ly, alpha, beta = struct.unpack("<ddfiiddff", raw_meta)
    meta = {"t": t, "dt": dt, "dtold": dtold, "nx": nx, "ny": ny,
            "Lx": lx, "Ly": ly, "alpha_vis": alpha, "beta_vis": beta}
    records = np.memmap(path, mode="r", dtype=record_dtype(), offset=16, shape=(count,))
    ids = records["u4if"] & np.uint64((1 << 56) - 1)
    if not np.array_equal(np.sort(ids), np.arange(count, dtype=np.uint64)):
        raise ValueError(f"Invalid or duplicate particle IDs in {path}")
    return records, ids, meta


def compare(reference_path, candidate_path):
    ref, ref_ids, ref_meta = load_snapshot(reference_path)
    cand, cand_ids, cand_meta = load_snapshot(candidate_path)
    if len(ref) != len(cand) or ref_meta["nx"] != cand_meta["nx"] or ref_meta["ny"] != cand_meta["ny"]:
        raise ValueError("Snapshot dimensions differ")
    ref_order = np.argsort(ref_ids)
    cand_order = np.argsort(cand_ids)
    if not np.array_equal(ref_ids[ref_order], cand_ids[cand_order]):
        raise ValueError("Particle ID sets differ")

    result = {
        "reference": str(Path(reference_path).resolve()),
        "candidate": str(Path(candidate_path).resolve()),
        "cells": len(ref),
        "record_size": RECORD_SIZE,
        "metadata": {"reference": ref_meta, "candidate": cand_meta,
                     "t_equal": ref_meta["t"] == cand_meta["t"]},
        "ids_equal": True,
        "type_tags_equal": bool(np.array_equal(ref["indt"][ref_order], cand["indt"][cand_order])),
        "primary_state_exact": True,
        "primary_state": {},
        "all_numeric_finite": True,
        "numeric_differences": {},
        "serialized_pointer_differences": {},
    }
    for name in ("next", "bp"):
        result["serialized_pointer_differences"][name] = int(
            np.count_nonzero(ref[name][ref_order] != cand[name][cand_order])
        )

    for field in PRIMARY_STATE:
        a = ref[field][ref_order]
        b = cand[field][cand_order]
        exact = bool(np.array_equal(a, b, equal_nan=True))
        result["primary_state"][field] = {"exact": exact}
        result["primary_state_exact"] &= exact

    ignored = {"next", "bp", "u4if", "indt"}
    for field in record_dtype().names:
        if field in ignored:
            continue
        a = ref[field][ref_order].astype(np.float64)
        b = cand[field][cand_order].astype(np.float64)
        finite = np.isfinite(a) & np.isfinite(b)
        result["all_numeric_finite"] &= bool(np.all(finite))
        if np.any(finite):
            delta = np.abs(a[finite] - b[finite])
            scale = np.maximum(np.abs(a[finite]), np.abs(b[finite]))
            field_scale = float(np.max(scale, initial=0.0))
            max_abs = float(np.max(delta, initial=0.0))
            rel_to_field_max = max_abs / field_scale if field_scale else 0.0
            rms = float(np.sqrt(np.mean(delta * delta)))
            changed = int(np.count_nonzero(delta))
        else:
            field_scale = 0.0
            max_abs = rel_to_field_max = rms = 0.0
            changed = 0
        if changed:
            result["numeric_differences"][field] = {
                "changed_cells": changed,
                "max_abs": max_abs,
                "rms_abs": rms,
                "max_field_abs": field_scale,
                "max_abs_over_field_max": rel_to_field_max,
            }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = compare(args.reference, args.candidate)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        with args.output.open("x") as stream:
            stream.write(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()

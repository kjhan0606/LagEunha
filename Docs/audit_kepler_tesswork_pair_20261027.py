"""Compare the accepted and instrumented one-step Kepler restarts by site ID."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("instrumented", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    p0, m0 = ordered(args.checkpoint, 55931)
    a, ma = ordered(args.baseline, 55931)
    b, mb = ordered(args.instrumented, 55931)
    if not (m0["t"] == 2.55 and ma["t"] == mb["t"] == 2.551 and
            np.array_equal(p0["mass"], a["mass"]) and
            np.array_equal(p0["mass"], b["mass"])):
        raise ValueError("Mismatched checkpoint, endpoint time, or fixed masses")

    fields = {}
    for name in a.dtype.names:
        xa, xb = a[name], b[name]
        unequal = xa != xb
        row = {"unequal_sites": int(np.count_nonzero(unequal)),
               "exact": bool(not np.any(unequal))}
        if np.issubdtype(xa.dtype, np.floating):
            row["maximum_absolute_difference"] = float(
                np.max(np.abs(xa.astype(float) - xb.astype(float))))
        fields[name] = row

    raw_a = np.memmap(args.baseline, mode="r", dtype=np.uint8)
    raw_b = np.memmap(args.instrumented, mode="r", dtype=np.uint8)
    if raw_a.size != raw_b.size:
        raise ValueError("Different snapshot sizes")
    diff = np.flatnonzero(raw_a != raw_b)
    record_start, record_end = 16, 16 + 55931 * 576
    record_diff = diff[(diff >= record_start) & (diff < record_end)]
    positions, counts = np.unique((record_diff - record_start) % 576,
                                   return_counts=True)
    order = np.argsort(counts)[::-1][:12]
    raw = {
        "sizes_equal": True,
        "unequal_bytes": int(diff.size),
        "header_unequal_bytes": int(np.sum(diff < record_start)),
        "record_unequal_bytes": int(record_diff.size),
        "footer_unequal_bytes": int(np.sum(diff >= record_end)),
        "first_unequal_byte": int(diff[0]) if diff.size else None,
        "most_changed_record_byte_offsets": [
            {"offset": int(positions[k]), "cells": int(counts[k])}
            for k in order
        ],
    }
    worst_ids = []
    for site_id in (46912, 55899):
        x, y = float(p0["x"][site_id]), float(p0["y"][site_id])
        worst_ids.append({
            "site_id": site_id, "x_at_t255": x, "y_at_t255": y,
            "radius_at_t255": float(np.hypot(x-3., y-3.)),
            "box_edge_distance_at_t255": min(x, y, 6.-x, 6.-y),
            "mass": float(p0["mass"][site_id]),
            "volume": float(p0["volume"][site_id]),
        })
    result = {
        "checkpoint": str(args.checkpoint.resolve()),
        "checkpoint_sha256": digest(args.checkpoint),
        "baseline": str(args.baseline.resolve()),
        "baseline_sha256": digest(args.baseline),
        "instrumented": str(args.instrumented.resolve()),
        "instrumented_sha256": digest(args.instrumented),
        "time": 2.551,
        "parsed_fields_all_exact": all(v["exact"] for v in fields.values()),
        "metadata_exact": ma == mb,
        "fields": fields,
        "raw_byte_comparison": raw,
        "kernel_work_extreme_site_positions": worst_ids,
        "scope": "One accepted step only; parsed fields omit portions of the 576-byte record.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in (
        "parsed_fields_all_exact", "metadata_exact", "raw_byte_comparison",
        "kernel_work_extreme_site_positions")}, indent=2))


if __name__ == "__main__":
    main()

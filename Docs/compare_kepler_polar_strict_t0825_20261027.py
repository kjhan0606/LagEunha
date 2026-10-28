"""Compare independently audited q1 and q2 native Kepler states at t=.825."""

import argparse
import hashlib
import json
from pathlib import Path


SCALARS = (
    "disk_density_abs_mass_weighted",
    "disk_entropy_log_abs_mass_weighted",
    "disk_entropy_log_signed_mass_weighted",
    "disk_position_error_abs_mass_weighted",
    "disk_angular_momentum_relative",
    "gizmo_L1_velocity",
    "gizmo_L1_density_literal",
    "gizmo_L1_density_plateau_corrected",
    "total_energy_relative",
)
BIN_METRICS = (
    "density_abs_mass_weighted",
    "entropy_log_abs_mass_weighted",
    "radial_speed_abs_mass_weighted",
    "position_error_abs_mass_weighted",
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def paired(a, b, signed=False):
    x, y = float(a), float(b)
    return {"q1": x, "q2": y,
            "q2_minus_q1": y - x,
            "q2_over_q1": None if signed or x == 0 else y / x}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("q1", type=Path)
    ap.add_argument("q2", type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    a = json.loads(args.q1.read_text())
    b = json.loads(args.q2.read_text())
    if (a["n"] != 384*384 or b["n"] != 384*384 or
            abs(a["time"] - .825) > 1e-12 or
            abs(b["time"] - .825) > 1e-12 or
            not all(a["checks"].values()) or not all(b["checks"].values())):
        raise ValueError("Require passing native q1/q2 strict receipts at exact t=.825")
    xbins = {(row["rlo"], row["rhi"]): row for row in a["bins"]}
    ybins = {(row["rlo"], row["rhi"]): row for row in b["bins"]}
    if xbins.keys() != ybins.keys():
        raise ValueError("Native q1/q2 radial bins differ")
    bins = []
    for key in sorted(xbins):
        x, y = xbins[key], ybins[key]
        bins.append({"rlo": key[0], "rhi": key[1],
                     "cells": {"q1": x["cells"], "q2": y["cells"]},
                     "metrics": {name: paired(x[name], y[name], "signed" in name)
                                 for name in BIN_METRICS}})
    result = {
        "time": .825,
        "q1_receipt": str(args.q1.resolve()), "q1_receipt_sha256": digest(args.q1),
        "q2_receipt": str(args.q2.resolve()), "q2_receipt_sha256": digest(args.q2),
        "n": 384*384,
        "scalars": {name: paired(a[name], b[name],
                                 "signed" in name or "relative" in name)
                    for name in SCALARS},
        "radial_bins": bins,
        "scope": "Matched-time native q1/q2 IC comparison; ratios do not isolate a single face or heating mechanism",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

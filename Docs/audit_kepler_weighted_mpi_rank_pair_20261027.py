"""Compare same-binary weighted runs at one versus two MPI ranks."""

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest

N = 55351
FIELDS = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2")


def runtime(path):
    text = path.read_text(errors="replace")
    match = re.search(r"\[RUN_TIME\]\s+wall=([0-9.eE+-]+)", text)
    if not match:
        raise ValueError(f"Missing native runtime marker: {path}")
    value = float(match.group(1))
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"Invalid runtime in {path}")
    return value


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("control", type=Path, help="one-rank audit receipt")
    ap.add_argument("candidate", type=Path, help="two-rank audit receipt")
    ap.add_argument("--binary-sha256", required=True)
    ap.add_argument("--parent-time", type=float, required=True)
    ap.add_argument("--time", type=float, required=True)
    ap.add_argument("--control-runtime", type=Path, required=True)
    ap.add_argument("--candidate-runtime", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    receipts = {"one_rank": json.loads(args.control.read_text()),
                "two_rank": json.loads(args.candidate.read_text())}
    for label, receipt in receipts.items():
        if (len(receipt["checks"]) != 21 or not all(receipt["checks"].values())
                or receipt["binary_sha256"] != args.binary_sha256
                or abs(receipt["parent_time"] - args.parent_time) > 1e-12
                or abs(receipt["time"] - args.time) > 1e-12):
            raise ValueError(f"{label} native physics, binary or time gate failed")
    a, b = receipts["one_rank"], receipts["two_rank"]
    if (a["initial_sha256"] != b["initial_sha256"]
            or a["parent_sha256"] != b["parent_sha256"]):
        raise ValueError("MPI runs do not share the same initial and parent states")

    arrays = {}
    for label, receipt in receipts.items():
        for phase in ("initial", "parent", "endpoint"):
            path = Path(receipt[phase])
            if digest(path) != receipt[phase + "_sha256"]:
                raise ValueError(f"Changed {label} {phase} snapshot")
            arrays[label, phase], _ = ordered(path, N)

    comparisons = {}
    for phase in ("initial", "parent", "endpoint"):
        comparisons[phase] = {}
        for field in FIELDS:
            x = arrays["one_rank", phase][field].astype(np.float64)
            y = arrays["two_rank", phase][field].astype(np.float64)
            delta = y - x
            rms = float(np.sqrt(np.mean(delta * delta)))
            scale = max(float(np.sqrt(np.mean(x * x))), float(np.max(np.abs(x))), 1e-300)
            comparisons[phase][field] = {
                "bitwise_equal": bool(np.array_equal(arrays["one_rank", phase][field],
                                                     arrays["two_rank", phase][field])),
                "max_abs_difference": float(np.max(np.abs(delta))),
                "rms_abs_difference": rms,
                "rms_relative_to_field_scale": rms / scale,
            }

    one_wall = runtime(args.control_runtime)
    two_wall = runtime(args.candidate_runtime)
    checks = {
        "both_native_physics_gates": True,
        "same_initial_state": True,
        "same_parent_state": True,
        "all_endpoint_fields_compared": all(field in comparisons["endpoint"] for field in FIELDS),
    }
    result = {
        "one_rank_receipt": str(args.control.resolve()),
        "two_rank_receipt": str(args.candidate.resolve()),
        "one_rank_receipt_sha256": digest(args.control),
        "two_rank_receipt_sha256": digest(args.candidate),
        "binary_sha256": args.binary_sha256,
        "parent_time": args.parent_time,
        "time": args.time,
        "ranks": {"control": 1, "candidate": 2},
        "runtime_seconds": {"one_rank": one_wall, "two_rank": two_wall,
                             "speedup_two_over_one": one_wall / two_wall},
        "fields": comparisons,
        "checks": checks,
        "scope": "Short same-checkpoint MPI scaling probe. Independent native gates are required; "
                 "rank-dependent floating-point differences are measured, not required to be bitwise zero.",
    }
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
    print(json.dumps({"checks": checks, "runtime_seconds": result["runtime_seconds"],
                      "endpoint_max_abs": {k: v["max_abs_difference"]
                                           for k, v in comparisons["endpoint"].items()}}, indent=2))


if __name__ == "__main__":
    main()

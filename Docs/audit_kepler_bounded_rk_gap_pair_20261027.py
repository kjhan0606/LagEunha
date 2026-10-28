"""Check paired native physics and exclusive RK4 timing intervals."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np


FIELDS = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2")
PROFILE_KEYS = ("wall", "tracked", "prelude", "after_update", "after_force",
                "after_post", "gaps", "closure")


def profile_rows(path):
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("[RK4 gap] "):
            continue
        values = {}
        for token in line.split()[2:]:
            key, value = token.split("=", 1)
            values[key] = int(value) if key == "step" else float(value.removesuffix("s"))
        if any(key not in values for key in ("step", *PROFILE_KEYS)):
            raise ValueError(f"Incomplete RK4 timing row: {line}")
        rows.append(values)
    if len(rows) < 3:
        raise ValueError("Need at least three RK4 profile rows")
    if len({row["step"] for row in rows}) != len(rows):
        raise ValueError("Repeated RK4 step number")
    for row in rows:
        if any(row[key] < 0 for key in PROFILE_KEYS[:-1]):
            raise ValueError("Negative measured interval")
        if abs(row["closure"]) >= .01:
            raise ValueError(f"Timing intervals do not close: {row}")
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("control", type=Path)
    parser.add_argument("profile", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    control = args.control.resolve()
    profile = args.profile.resolve()
    sys.path.insert(0, str(control))
    from audit_kepler_sparse_area_t0825_20261027 import ordered
    from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time

    receipts = [json.loads((directory / "audit_t085.json").read_text())
                for directory in (control, profile)]
    if not all(all(receipt["checks"].values()) for receipt in receipts):
        raise ValueError("A paired native physics gate failed")
    if receipts[0]["restart_sha256"] != receipts[1]["restart_sha256"]:
        raise ValueError("Different accepted restart snapshots")
    starts = [directory / "khout.000001.dat" for directory in (control, profile)]
    endpoints = [Path(receipt["endpoint"]) for receipt in receipts]
    timed_paths = list(zip(starts, (.825, .825))) + list(zip(endpoints, (.85, .85)))
    if any(abs(snapshot_time(path, 55351) - target) > 1e-10
           for path, target in timed_paths):
        raise ValueError("Unexpected paired checkpoint time")
    states = [ordered(path, 55351)[0] for path in endpoints]
    differences = {}
    for key in FIELDS:
        left, right = states[0][key], states[1][key]
        differences[key] = {
            "bitwise_equal": bool(np.array_equal(left, right)),
            "absolute_max": float(np.max(np.abs(left.astype(float) - right.astype(float)))),
        }
    rows = profile_rows(profile / "log")
    steady = rows[1:]
    medians = {key: float(np.median([row[key] for row in steady]))
               for key in PROFILE_KEYS}
    result = {
        "control": str(control), "profile": str(profile),
        "restart_sha256": receipts[0]["restart_sha256"],
        "control_endpoint_sha256": digest(endpoints[0]),
        "profile_endpoint_sha256": digest(endpoints[1]),
        "physical_fields": differences,
        "profile_rows": rows,
        "median_excluding_first_step_seconds": medians,
        "scope": "Paired short native restart; timing bins are exclusive within RK4, GPU tessellation is inside update",
    }
    rendered = json.dumps(result, indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as stream:
        stream.write(rendered)
    print(rendered)
    if not all(row["bitwise_equal"] for row in differences.values()):
        raise SystemExit("Profile patch changed a parsed physical field")


if __name__ == "__main__":
    main()

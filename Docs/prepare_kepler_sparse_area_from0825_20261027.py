"""Promote corrected sparse Kepler t=.825 only after every physical gate."""

import argparse
import json
from pathlib import Path
import re
import shutil

from prepare_kepler_sparse_area_restart_20261027 import (
    BIN_SHA, IC, digest, snapshot_time,
)


TARGET = 10.
REQUIRED = (
    "initial_periodic_volume", "endpoint_periodic_volume", "mass_exact",
    "W0_exact", "positive_state", "no_floor_injection",
    "total_energy_relative_below_1e_minus_3", "disk_density_abs_below_0p005",
    "disk_entropy_abs_below_0p01", "disk_position_error_below_1e_minus_4",
    "disk_radial_speed_rms_below_0p01",
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("predecessor", type=Path)
    ap.add_argument("destination", type=Path)
    ap.add_argument("--expected-count", type=int, choices=tuple(IC), required=True)
    args = ap.parse_args()
    parent = args.predecessor.resolve()
    dest = args.destination.resolve()
    receipt = json.loads((parent / "audit_t0825.json").read_text())
    if (receipt["n"] != args.expected_count or abs(receipt["time"] - .825) > 1e-12 or
            not all(receipt["checks"].values()) or
            not all(receipt["checks"].get(key) is True for key in REQUIRED) or
            abs(receipt["initial_volume_sum"] - 36) > 1e-6 or
            abs(receipt["endpoint_volume_sum"] - 36) > 1e-5):
        raise ValueError("Corrected sparse t=.825 physical audit did not pass")
    initial = Path(receipt["initial"])
    source = Path(receipt["endpoint"])
    if source.parent.resolve() != parent or not initial.is_file():
        raise ValueError("Audit checkpoint lineage is invalid")
    preflight = initial.parent.resolve()
    first = json.loads((preflight / "audit_t001.json").read_text())
    if (first["n"] != args.expected_count or
            digest(initial) != first["sha256"].get(initial.name) or
            abs(first["initial_total_mass"] / receipt["mass_sum"] - 1) > 1e-10):
        raise ValueError("Corrected preflight initial state changed")
    match = re.fullmatch(r"khout\.(\d+)\.dat", source.name)
    if (not match or int(match.group(1)) <= 0 or
            abs(snapshot_time(source, args.expected_count) - .825) > 1e-12 or
            digest(source) != receipt["endpoint_sha256"] or
            digest(parent / "eunha2.exe") != BIN_SHA):
        raise ValueError("t=.825 checkpoint or accepted binary changed")
    if dest.exists() and (any(dest.glob("khout.*.dat")) or
                          any((dest / name).exists() for name in
                              ("eunha2.exe", "params.dat", "restart_plan.json"))):
        raise ValueError("Destination already has native state inputs")
    dest.mkdir(parents=True, exist_ok=True)
    copied = {}
    for original in (source, parent / "eunha2.exe", parent / "params.dat"):
        target = dest / original.name
        shutil.copy2(original, target)
        copied[target.name] = digest(target)
    plan = {"predecessor_receipt": str(parent / "audit_t0825.json"),
            "predecessor_checkpoint": str(source),
            "predecessor_sha256": receipt["endpoint_sha256"],
            "initial": str(initial), "step": int(match.group(1)),
            "expected_time": .825, "target_time": TARGET,
            "n": args.expected_count, "binary_sha256": BIN_SHA,
            "inputs_sha256": copied}
    with (dest / "restart_plan.json").open("x") as stream:
        json.dump(plan, stream, indent=2)
        stream.write("\n")
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()

"""Restart q1 from its last complete audited snapshot after halo clipping."""

import argparse
import json
from pathlib import Path
import re
import shutil

from prepare_kepler_polar_long_restart_20261027 import (
    BIN_SHA256, N, digest, require_zero_injection, snapshot_time)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("parent", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("initial", type=Path)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    if abs(snapshot_time(args.initial)) > 1e-12:
        raise ValueError("Initial q1 reference is not exact t=0")
    log = (args.parent / "log").read_text(errors="replace")
    if ("[GPU_TESS_INVALID] volume=0 velocity_fit=0 geometry=0 halo=4" not in log
            or "[RK4_GEOMETRY_ABORT_CHECKPOINT]" not in log
            or "EXIT:1" not in log):
        raise ValueError("Parent is not the documented q1 halo-clip failure")
    require_zero_injection(args.parent / "log")
    receipt = json.loads(args.receipt.read_text())
    if (receipt["n"] != N or not all(receipt["checks"].values())
            or abs(receipt["time"] - .6187644470884642) > 1e-12):
        raise ValueError("Last complete q1 checkpoint audit did not pass")
    source = Path(receipt["endpoint"])
    match = re.fullmatch(r"khout\.(\d+)\.dat", source.name)
    if (not match or int(match.group(1)) != 8119
            or source.parent.resolve() != args.parent.resolve()
            or abs(snapshot_time(source) - receipt["time"]) > 1e-12):
        raise ValueError("Receipt does not identify complete accepted step 8119")
    if digest(args.parent / "eunha2.exe") != BIN_SHA256:
        raise ValueError("Parent q1 binary changed")
    if args.destination.exists() and any(args.destination.glob("khout.*.dat")):
        raise ValueError("Destination already contains evolution snapshots")
    args.destination.mkdir(parents=True, exist_ok=True)
    inputs = {}
    for original in (source, args.parent / "eunha2.exe", args.parent / "params.dat"):
        target = args.destination / original.name
        if target.exists():
            raise ValueError(f"Destination input already exists: {target}")
        shutil.copy2(original, target)
        inputs[target.name] = digest(target)
    result = {"parent_failure_log": str(args.parent / "log"),
              "accepted_audit": str(args.receipt),
              "source": str(source),
              "source_time": receipt["time"],
              "source_step": 8119,
              "target_time": .825,
              "new_halo_cells": 12,
              "inputs_sha256": inputs,
              "scope": "Independent q1 halo-width recovery experiment; never uses failed-stage khdiag as restart"}
    with (args.destination / "recovery_plan.json").open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

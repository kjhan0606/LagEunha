"""Copy one audited q1 or q2 native checkpoint into a new long-time restart run."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import struct


N = 384 * 384
STRIDE = 576
BIN_SHA256 = "7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def snapshot_time(path):
    with path.open("rb") as stream:
        magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
        stream.seek(-52, 2)
        time = struct.unpack("<d", stream.read(8))[0]
    if ((magic, n, stride, stress) != (-1, N, STRIDE, 1)
            or path.stat().st_size != 16 + N * STRIDE + 52):
        raise ValueError(f"Incomplete or incompatible checkpoint: {path}")
    return time


def require_zero_injection(log_path):
    lines = [line for line in log_path.read_text(errors="replace").splitlines()
             if "[RK4E]" in line]
    if not lines:
        raise ValueError("Parent native RK4 energy ledger is missing")
    last = lines[-1]
    for key in ("floor_cum", "de_cum", "sfl_cum"):
        match = re.search(rf"\b{key}=([0-9.eE+-]+)", last)
        if not match or float(match.group(1)) != 0.:
            raise ValueError(f"Parent {key} is missing or nonzero")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("parent", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("initial", type=Path)
    parser.add_argument("--receipt", required=True)
    parser.add_argument("--expected-time", required=True, type=float)
    parser.add_argument("--target-time", required=True, type=float)
    args = parser.parse_args()
    if not args.target_time > args.expected_time > 0:
        raise ValueError("Invalid restart time interval")
    if abs(snapshot_time(args.initial)) > 1e-12:
        raise ValueError("Initial reference is not exact t=0")
    receipt_path = args.parent / args.receipt
    receipt = json.loads(receipt_path.read_text())
    if (not all(receipt["checks"].values()) or receipt["n"] != N
            or abs(receipt["time"] - args.expected_time) > 1e-10):
        raise ValueError("Parent native physics gate did not pass")
    require_zero_injection(args.parent / "log")
    source = Path(receipt["endpoint"])
    match = re.fullmatch(r"khout\.(\d+)\.dat", source.name)
    if (not match or int(match.group(1)) <= 0
            or source.parent.resolve() != args.parent.resolve()
            or abs(snapshot_time(source) - args.expected_time) > 1e-10):
        raise ValueError("Receipt checkpoint is not the complete expected parent")
    if digest(args.parent / "eunha2.exe") != BIN_SHA256:
        raise ValueError("Parent executable does not match accepted q1/q2 binary")
    if args.destination.exists() and any(args.destination.glob("khout.*.dat")):
        raise ValueError("Destination already has an evolution checkpoint")
    args.destination.mkdir(parents=True, exist_ok=True)
    inputs = {}
    for original in (source, args.parent / "eunha2.exe", args.parent / "params.dat"):
        target = args.destination / original.name
        if target.exists():
            raise ValueError(f"Destination input already exists: {target}")
        shutil.copy2(original, target)
        inputs[target.name] = digest(target)
    plan = {"predecessor_receipt": str(receipt_path),
            "predecessor_checkpoint": str(source),
            "expected_time": args.expected_time,
            "target_time": args.target_time,
            "step": int(match.group(1)),
            "initial": str(args.initial),
            "inputs_sha256": inputs}
    with (args.destination / "restart_plan.json").open("x") as stream:
        json.dump(plan, stream, indent=2)
        stream.write("\n")
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()

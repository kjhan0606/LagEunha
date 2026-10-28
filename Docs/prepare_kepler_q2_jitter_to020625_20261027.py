"""Promote only a passing site-jitter preflight to the q2 comparison time."""

import argparse
import json
from pathlib import Path
import re
import shutil

from prepare_kepler_polar_long_restart_20261027 import (
    BIN_SHA256, N, digest, require_zero_injection, snapshot_time,
)


JITTER_SHA256 = "754c6c60acb516404f6799111dac6605ad58c04baa166f706144c046e6a5d080"
# The accepted q2 snapshot khout.002700.dat is at this actual saved time,
# slightly beyond its nominal 0.20625 dump clock.
TARGET = 0.20626460064652694


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("preflight", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    parent = args.preflight.resolve()
    destination = args.destination.resolve()
    audit = json.loads((parent / "audit_t001.json").read_text())
    if (audit["n"] != N or abs(audit["t0"]) > 1e-12
            or abs(audit["t1"] - .001) > 1e-12
            or not all(audit["checks"].values())):
        raise ValueError("Jitter native preflight did not pass")
    if (digest(parent / "glass.bin") != JITTER_SHA256 or
            digest(parent / "eunha2.exe") != BIN_SHA256):
        raise ValueError("Jitter glass or accepted native binary changed")
    initial = Path(audit["initial"])
    source = Path(audit["evolved"])
    if initial.parent.resolve() != parent or source.parent.resolve() != parent:
        raise ValueError("Preflight receipt points outside its run")
    if abs(snapshot_time(initial)) > 1e-12 or abs(snapshot_time(source) - .001) > 1e-12:
        raise ValueError("Preflight checkpoint time/header is wrong")
    match = re.fullmatch(r"khout\.(\d+)\.dat", source.name)
    if not match or int(match.group(1)) <= 0:
        raise ValueError("Unexpected preflight checkpoint name")
    for path in (initial, source):
        if digest(path) != audit["sha256"][path.name]:
            raise ValueError("Preflight checkpoint changed after its audit")
    require_zero_injection(parent / "log")
    if destination.exists() and (
            any(destination.glob("khout.*.dat")) or
            any((destination / name).exists() for name in
                ("eunha2.exe", "params.dat", "restart_plan.json"))):
        raise ValueError("Restart destination already has state inputs")
    destination.mkdir(parents=True, exist_ok=True)
    copied = {}
    for original in (source, parent / "eunha2.exe", parent / "params.dat"):
        target = destination / original.name
        shutil.copy2(original, target)
        copied[target.name] = digest(target)
    plan = {"preflight_receipt": str(parent / "audit_t001.json"),
            "preflight_checkpoint": str(source),
            "initial": str(initial), "step": int(match.group(1)),
            "expected_time": .001, "target_time": TARGET,
            "jitter_glass_sha256": JITTER_SHA256,
            "binary_sha256": BIN_SHA256,
            "inputs_sha256": copied}
    with (destination / "restart_plan.json").open("x") as stream:
        json.dump(plan, stream, indent=2)
        stream.write("\n")
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()

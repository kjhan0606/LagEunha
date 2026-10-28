"""Copy a corrected sparse Kepler restart only from audited native state."""

import argparse
import json
import math
from pathlib import Path
import re
import shutil

from prepare_kepler_sparse_area_restart_20261027 import (
    BIN_SHA, IC, digest, snapshot_time,
)


PARAMS_SHA = "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("predecessor", type=Path)
    ap.add_argument("destination", type=Path)
    ap.add_argument("--receipt", required=True)
    ap.add_argument("--consistency-receipt", required=True)
    ap.add_argument("--expected-count", type=int, choices=tuple(IC), required=True)
    ap.add_argument("--expected-time", type=float, required=True)
    ap.add_argument("--target-time", type=float, required=True)
    args = ap.parse_args()
    parent = args.predecessor.resolve()
    dest = args.destination.resolve()
    if (dest == parent or not math.isfinite(args.expected_time) or
            not math.isfinite(args.target_time) or args.expected_time < 10. or
            args.target_time <= args.expected_time):
        raise ValueError("Expected later target and an independent destination")
    audit_path = parent / args.receipt
    consistency_path = parent / args.consistency_receipt
    audit = json.loads(audit_path.read_text())
    consistency = json.loads(consistency_path.read_text())
    n, t = args.expected_count, args.expected_time
    if (audit.get("n") != n or consistency.get("n") != n or
            abs(float(audit["time"]) - t) > 1e-10 or
            abs(float(consistency["time"]) - t) > 1e-10 or
            not all(audit["checks"].values()) or
            not all(consistency["checks"].values())):
        raise ValueError("Predecessor physics/state audits did not pass")
    initial = Path(audit["initial"]).resolve()
    source = Path(audit["endpoint"]).resolve()
    if (source.parent != parent or
            initial != Path(consistency["initial"]).resolve() or
            source != Path(consistency["endpoint"]).resolve()):
        raise ValueError("Audit checkpoint lineage disagrees")
    match = re.fullmatch(r"khout\.(\d+)\.dat", source.name)
    if (not match or int(match.group(1)) <= 0 or
            abs(snapshot_time(initial, n)) > 1e-12 or
            abs(snapshot_time(source, n) - t) > 1e-10):
        raise ValueError("Initial or restart snapshot is incomplete or mistimed")
    initial_sha = digest(initial)
    source_sha = digest(source)
    if (initial_sha != consistency["initial_sha256"] or
            source_sha != consistency["endpoint_sha256"] or
            source_sha != audit["endpoint_sha256"] or
            digest(parent / "eunha2.exe") != BIN_SHA or
            digest(parent / "params.dat") != PARAMS_SHA):
        raise ValueError("Audited checkpoint, binary, or physical input changed")
    if dest.exists() and (any(dest.glob("khout.*.dat")) or
                          any((dest / name).exists() for name in
                              ("eunha2.exe", "params.dat", "restart_plan.json"))):
        raise ValueError("Destination already contains a native restart")
    dest.mkdir(parents=True, exist_ok=True)
    copied = {}
    for original in (source, parent / "eunha2.exe", parent / "params.dat"):
        target = dest / original.name
        shutil.copy2(original, target)
        copied[target.name] = digest(target)
    plan = {
        "predecessor_receipt": str(audit_path),
        "predecessor_consistency_receipt": str(consistency_path),
        "predecessor_checkpoint": str(source),
        "predecessor_checkpoint_sha256": source_sha,
        "initial": str(initial), "initial_sha256": initial_sha,
        "step": int(match.group(1)), "n": n,
        "expected_time": t, "target_time": args.target_time,
        "binary_sha256": BIN_SHA, "params_sha256": PARAMS_SHA,
        "inputs_sha256": copied,
    }
    with (dest / "restart_plan.json").open("x") as stream:
        json.dump(plan, stream, indent=2)
        stream.write("\n")
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()

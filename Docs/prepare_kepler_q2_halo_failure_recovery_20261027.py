"""After a q2 halo abort, select the latest verified accepted RK4 state."""

import argparse
import json
import math
from pathlib import Path
import re
import shutil
import struct
import sys

import numpy as np

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot

from prepare_kepler_polar_long_restart_20261027 import (
    BIN_SHA256, N, digest, require_zero_injection, snapshot_time)


FIELDS = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2")
DIAG_RE = (r"\[RK4_FAILURE_SNAPSHOT\] "
           r"(khdiag\.step(\d+)\.start_geometry\.dat) np=(\d+) "
           r"t_start=([0-9.eE+-]+); failed stage is not a restart state")


def ordered_if_valid(path):
    try:
        t = snapshot_time(path)
        p, meta = read_snapshot(path)
        ids = (p["u4if"].astype(np.uint64)
               & np.uint64((1 << 56) - 1)).astype(np.int64)
        if (len(p) != N or not np.array_equal(np.sort(ids), np.arange(N))
                or abs(float(meta["t"]) - t) > 1e-12
                or not all(np.isfinite(p[key]).all() for key in FIELDS)):
            return None
        return p[np.argsort(ids)], t
    except (OSError, ValueError, KeyError, struct.error):
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("parent", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("initial", type=Path)
    args = parser.parse_args()
    log = (args.parent / "log").read_text(errors="replace")
    failure = re.search(r"\[GPU_TESS_INVALID\] volume=0 velocity_fit=0 geometry=0 halo=([1-9][0-9]*)", log)
    if not failure or "[RK4_GEOMETRY_ABORT_CHECKPOINT]" not in log or "EXIT:1" not in log:
        raise ValueError("Parent did not end with the expected isolated halo abort")
    require_zero_injection(args.parent / "log")
    initial = ordered_if_valid(args.initial)
    if initial is None or initial[1] != 0:
        raise ValueError("Initial q2 snapshot failed completeness check")
    p0, _ = initial
    if (not np.all(p0["w2"] == 0) or not np.all(p0["ie"] > 0)
            or not np.all(p0["volume"] > 0)):
        raise ValueError("Initial q2 physical state invalid")
    if digest(args.parent / "eunha2.exe") != BIN_SHA256:
        raise ValueError("Parent q2 binary changed")
    candidates = []
    for path in args.parent.glob("khout.*.dat"):
        match = re.fullmatch(r"khout\.(\d+)\.dat", path.name)
        state = ordered_if_valid(path)
        if not match or state is None:
            continue
        p, t = state
        if not (.001 - 1e-12 <= t < .825 - 1e-12):
            continue
        if (not np.array_equal(p["mass"], p0["mass"])
                or not np.all(p["w2"] == 0)
                or not np.all(p["ie"] > 0)
                or not np.all(p["volume"] > 0)):
            continue
        candidates.append((t, int(match.group(1)), path))
    if not candidates:
        raise ValueError("No complete positive fixed-mass W=0 q2 checkpoint")
    ordinary_time, ordinary_step, ordinary_source = max(candidates)
    t, step, source = ordinary_time, ordinary_step, ordinary_source
    source_kind = "ordinary_accepted_output"
    diagnostic = re.search(DIAG_RE, log)
    if diagnostic:
        name, diag_step, diag_count, logged_time = diagnostic.groups()
        diag_step, diag_count, logged_time = int(diag_step), int(diag_count), float(logged_time)
        if (diag_count != N or diag_step <= ordinary_step or
                "[RK4_GEOMETRY_ABORT_CHECKPOINT] accepted step-start saved; failed geometry is not serialized" not in log):
            raise ValueError("Geometry-abort state is not a verified accepted step start")
        candidate = args.parent / name
        state = ordered_if_valid(candidate)
        if state is None:
            raise ValueError("Accepted geometry-abort state is incomplete")
        pdiag, diag_time = state
        with candidate.open("rb") as stream:
            stream.seek(-52, 2)
            trailer = struct.unpack("<ddfiiddff", stream.read(52))
        metadata_ok = (abs(diag_time-logged_time) < 1e-12 and
                       abs(trailer[0]-diag_time) < 1e-12 and
                       all(math.isfinite(value) and value > 0 for value in trailer[1:3]) and
                       trailer[3:7] == (384, 384, 6., 6.) and
                       .001 < ordinary_time < diag_time < .825)
        physical_ok = (np.array_equal(pdiag["mass"], p0["mass"]) and
                       np.all(pdiag["w2"] == 0) and
                       np.all(pdiag["ie"] > 0) and np.all(pdiag["volume"] > 0) and
                       abs(float(np.sum(pdiag["volume"], dtype=np.longdouble))-36.) < 1e-5 and
                       np.max(abs(pdiag["den"].astype(float) /
                                  (pdiag["mass"].astype(float)/pdiag["volume"].astype(float))-1.)) < 1e-6)
        if not (metadata_ok and physical_ok):
            raise ValueError("Accepted geometry-abort state failed restart metadata or physics gate")
        t, step, source = diag_time, diag_step, candidate
        source_kind = "verified_accepted_k1_start"
    if source.parent.resolve() != args.parent.resolve():
        raise ValueError("Selected checkpoint escaped the failed run")
    if args.destination.exists() and (any(args.destination.glob("khout.*.dat")) or
                                      any((args.destination / name).exists() for name in
                                          ("eunha2.exe", "params.dat", "recovery_plan.json"))):
        raise ValueError("Destination already contains native restart inputs")
    args.destination.mkdir(parents=True, exist_ok=True)
    inputs = {}
    for original in (source, args.parent / "eunha2.exe", args.parent / "params.dat"):
        target = args.destination / (f"khout.{step:06d}.dat" if original == source else original.name)
        if target.exists():
            raise ValueError(f"Destination input already exists: {target}")
        shutil.copy2(original, target)
        inputs[target.name] = digest(target)
    plan = {"parent_failure_log": str(args.parent / "log"),
            "halo_clipped_cells": int(failure.group(1)),
            "source": str(source), "source_kind": source_kind,
            "source_sha256": digest(source),
            "last_ordinary_source": str(ordinary_source),
            "last_ordinary_time": ordinary_time,
            "source_time": t, "source_step": step,
            "target_time": .825, "new_halo_cells": 12,
            "initial": str(args.initial), "inputs_sha256": inputs,
            "scope": "Conditional q2 halo-width recovery from byte-audited accepted RK4 state; failed-stage geometry excluded"}
    with (args.destination / "recovery_plan.json").open("x") as stream:
        json.dump(plan, stream, indent=2)
        stream.write("\n")
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()

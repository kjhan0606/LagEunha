#!/usr/bin/env python3
"""Stage a one-step CUDA polygon-moment audit from an accepted q1 restart."""

from __future__ import annotations

import hashlib
import json
import shutil
import struct
import subprocess
from pathlib import Path


SOURCE_ROOT = Path("/tmp/Kepler_spin_moments_20261028")
BASE = Path("/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1")
PARENT = BASE / "dteta_006_t11_probe_20261028_a10"
DESTINATION = BASE / "polygon_moment_audit_20261028_a10"
RESTART_STEP = 13200
RESTART_TIME = 1.0000125425182276
TARGET_TIME = RESTART_TIME + 5.0e-5
N = 384 * 384
STRIDE = 576
EXPECTED_BASE_BINARY = "7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284"
EXPECTED_PARAMS = "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"
EXPECTED_RESTART = "24ed4e62d0160f26659a71f490360c42e196e6592c1152db2ea06dec9f01bea8"
EXPECTED_INITIAL = "edb1f9d30700a95093a4f442afd3461238dd4e93ad9b075302557c15dcd693c9"
JOB_SCRIPT = Path(__file__).with_name("kepler_polygon_moment_audit_20261028.slurm")


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def check_snapshot(path: Path, expected_time: float) -> None:
    if path.stat().st_size != 16 + N * STRIDE + 52:
        raise ValueError(f"Unexpected snapshot size: {path}")
    with path.open("rb") as stream:
        header = struct.unpack("<4i", stream.read(16))
        stream.seek(-52, 2)
        actual_time = struct.unpack("<d", stream.read(8))[0]
    if header != (-1, N, STRIDE, 1) or abs(actual_time - expected_time) > 1.0e-10:
        raise ValueError(f"Wrong snapshot header/time for {path}: {header}, {actual_time}")


def main() -> None:
    if DESTINATION.exists():
        raise FileExistsError(f"Refusing to reuse {DESTINATION}")
    if not SOURCE_ROOT.is_dir() or not (SOURCE_ROOT / ".git").exists():
        raise FileNotFoundError(f"Expected isolated source worktree: {SOURCE_ROOT}")
    if not JOB_SCRIPT.is_file():
        raise FileNotFoundError(JOB_SCRIPT)

    plan = json.loads((PARENT / "restart_plan.json").read_text())
    checkpoint = PARENT / plan["checkpoint"]
    initial = Path(plan["initial"])
    params = PARENT / "params.dat"
    baseline_binary = PARENT / "eunha2.exe"
    if (plan["restart_step"] != RESTART_STEP or
            abs(plan["restart_time"] - RESTART_TIME) > 1.0e-12 or
            digest(checkpoint) != EXPECTED_RESTART or
            digest(initial) != EXPECTED_INITIAL or
            digest(params) != EXPECTED_PARAMS or
            digest(baseline_binary) != EXPECTED_BASE_BINARY):
        raise ValueError("Accepted q1 restart lineage changed; refusing to stage")
    check_snapshot(checkpoint, RESTART_TIME)
    check_snapshot(initial, 0.0)

    binary = SOURCE_ROOT / "eunha2.exe"
    if not binary.is_file() or binary.stat().st_size == 0:
        raise FileNotFoundError(f"Candidate executable missing: {binary}")
    source_hashes = {
        name: digest(SOURCE_ROOT / name)
        for name in ("Exam/exam_gpu.cu", "Exam/exam_gpu.h", "Exam/exam_gpu_extract.c",
                     "Exam/gfs_polygon_moments_reference_20261028.h")
    }
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT, text=True
    ).strip()

    DESTINATION.mkdir(parents=False, exist_ok=False)
    files = {
        binary: DESTINATION / "eunha2.exe",
        params: DESTINATION / "params.dat",
        checkpoint: DESTINATION / checkpoint.name,
        initial: DESTINATION / initial.name,
        JOB_SCRIPT: DESTINATION / "job.slurm",
    }
    for source, target in files.items():
        shutil.copy2(source, target)
    (DESTINATION / "job.slurm").chmod(0o755)

    candidate = {
        "purpose": "one-step GPU polygon area/centroid/mass-Q audit; no spin evolution",
        "source_root": str(SOURCE_ROOT),
        "source_commit": commit,
        "source_hashes": source_hashes,
        "binary_sha256": digest(DESTINATION / "eunha2.exe"),
        "baseline_binary_sha256": EXPECTED_BASE_BINARY,
        "params_sha256": digest(DESTINATION / "params.dat"),
        "checkpoint": checkpoint.name,
        "checkpoint_sha256": digest(DESTINATION / checkpoint.name),
        "initial": initial.name,
        "initial_sha256": digest(DESTINATION / initial.name),
        "restart_step": RESTART_STEP,
        "restart_time": RESTART_TIME,
        "target_time": TARGET_TIME,
        "dt_eta": 0.06,
        "geometry": "GFS_GPU_RAW_FACES=1; every local polygon compared once to the host long-double reference",
        "requested_resources": "1 MPI rank, 2 OpenMP threads, 1 A10 GPU",
        "dependency": "afterany:438007; exclude syn02",
    }
    with (DESTINATION / "restart_plan.json").open("x") as stream:
        json.dump(candidate, stream, indent=2)
        stream.write("\n")
    (DESTINATION / "OWNER.txt").write_text(
        "Project: LagForce / LagEunha\n"
        "Purpose: one-step CUDA polygon-moment geometry audit from the accepted q1 restart.\n"
        "Source: isolated /tmp/Kepler_spin_moments_20261028 worktree; see restart_plan.json.\n"
        "DO_NOT_CANCEL until the audit result is reviewed.\n"
    )
    checked = ["eunha2.exe", "params.dat", checkpoint.name, initial.name, "restart_plan.json"]
    (DESTINATION / "inputs.sha256").write_text("".join(
        f"{digest(DESTINATION / name)}  {name}\n" for name in checked
    ))
    print(json.dumps({
        "destination": str(DESTINATION),
        "binary_sha256": candidate["binary_sha256"],
        "checkpoint_sha256": candidate["checkpoint_sha256"],
        "restart_time": RESTART_TIME,
        "target_time": TARGET_TIME,
        "source_commit": commit,
        "source_hashes": source_hashes,
    }, indent=2))


if __name__ == "__main__":
    main()

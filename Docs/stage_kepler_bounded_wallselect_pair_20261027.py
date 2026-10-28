"""Stage a same-compiler A10 control of selective fixed-wall preclipping."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import snapshot_time


DOCS = Path(__file__).resolve().parent
BASE = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/to0825_near32_v1")
SOURCES = {
    "control": Path("/tmp/Kepler_sparse_fixedpoly_profile_20261027"),
    "candidate": Path("/tmp/Kepler_sparse_fixedpoly_wallselect_20261027"),
}
GPU_SOURCES = {
    "control": "2194c4aa81ea0eee030c5b07955794d0b7dc8d25dbc32f5049d14c71add67865",
    "candidate": "9011dd6de77c2e938400e0fd7bce7dbbbedbd671eee028dabf066f643b498a6e",
}
CPU_SOURCE = "6b1d54c9de7ed7ce710293d92daac7e055434acf09a3674dd4026b9cde0c58f9"
CHECKPOINT_SHA = "134e4e76a475c6923f1106ffb944fd0d51a0df82e7d4c86ea2814057ce9a12f6"
HELPERS = (
    "audit_kepler_bounded_short_replay_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)
RUNNER = "kepler_q2_bounded_wallselect_pair_t085_a10_20261027.slurm"
AUDITOR = "audit_kepler_bounded_rk_gap_pair_20261027.py"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    source_snapshot = BASE / "khout.000403.dat"
    if digest(source_snapshot) != CHECKPOINT_SHA or abs(snapshot_time(source_snapshot, 55351) - .825) > 1e-10:
        raise ValueError("Accepted checkpoint changed")
    receipt = json.loads((BASE / "audit_t0825.json").read_text())
    if receipt["endpoint_sha256"] != CHECKPOINT_SHA or not all(receipt["checks"].values()):
        raise ValueError("Predecessor physics gate failed")
    for name, source in SOURCES.items():
        if digest(source / "Exam/exam.c") != CPU_SOURCE:
            raise ValueError(f"{name} CPU source changed")
        if digest(source / "Exam/exam_gpu.cu") != GPU_SOURCES[name]:
            raise ValueError(f"{name} GPU source changed")
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)
    dest.mkdir(parents=True)
    for name in (RUNNER, AUDITOR):
        shutil.copy2(DOCS / name, dest / name)
    hashes = {}
    for variant, source in SOURCES.items():
        directory = dest / variant
        directory.mkdir()
        inputs = {
            "eunha2.exe": source / "eunha2.exe",
            "params.dat": BASE / "params.dat",
            "khout.000001.dat": source_snapshot,
            **{name: DOCS / name for name in HELPERS},
        }
        hashes[variant] = {}
        for name, path in inputs.items():
            shutil.copy2(path, directory / name)
            hashes[variant][name] = digest(directory / name)
        with (directory / "inputs.sha256").open("x") as stream:
            for name, sha in sorted(hashes[variant].items()):
                stream.write(f"{sha}  {name}\n")
    manifest = {
        "destination": str(dest), "source": str(BASE),
        "checkpoint_time": .825, "target_time": .85,
        "checkpoint_sha256": CHECKPOINT_SHA,
        "gpu_source_sha256": GPU_SOURCES,
        "cpu_source_sha256": CPU_SOURCE,
        "runner_sha256": digest(dest / RUNNER),
        "auditor_sha256": digest(dest / AUDITOR),
        "inputs": hashes,
        "scope": "Same-compiler short CUDA pair; only selective 24-neighbor wall preclip differs",
    }
    (dest / "stage_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

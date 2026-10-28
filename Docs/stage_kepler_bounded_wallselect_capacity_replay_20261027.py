"""Replay the accepted t=.0456165 corner-capacity state with the GPU candidate."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
SOURCE = Path("/tmp/Kepler_sparse_fixedpoly_wallselect_20261027/eunha2.exe")
SOURCE_SHA = "87d60c1cb059e1710dd2749c84af5c98ba17fe077e6fe080b9c8b42d5fd40968"
CHECKPOINT = ROOT / "to0825_v2/khdiag.step000024.start_geometry.dat"
CHECKPOINT_SHA = "08a7500e5d72c965c771ad51fd90d87a3cd7872a439edfb146df7a566d18c76a"
DEST = ROOT / "wallselect_capacity_replay005_a10_v1"
RUNNER = "kepler_q2_bounded_capacity_replay_to005_a10_20261027.slurm"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    if digest(SOURCE) != SOURCE_SHA or digest(CHECKPOINT) != CHECKPOINT_SHA:
        raise ValueError("Candidate binary or accepted checkpoint changed")
    if abs(snapshot_time(CHECKPOINT, 55351) - .045616533734626756) > 1e-10:
        raise ValueError("Wrong accepted checkpoint time")
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)
    files = {
        "eunha2.exe": SOURCE,
        "params.dat": ROOT / "t001_v2/params.dat",
        "khout.000001.dat": CHECKPOINT,
        "job.slurm": DOCS / RUNNER,
        **{name: DOCS / name for name in (
            "audit_kepler_bounded_short_replay_20261027.py",
            "audit_kepler_sparse_area_t0825_20261027.py",
            "audit_kepler_sparse_state_consistency_20261027.py",
            "prepare_kepler_sparse_area_restart_20261027.py",
        )},
    }
    hashes = {name: digest(path) for name, path in files.items()}
    DEST.mkdir()
    for name, path in files.items():
        shutil.copy2(path, DEST / name)
        if digest(DEST / name) != hashes[name]:
            raise IOError(f"Staged {name} changed")
    with (DEST / "inputs.sha256").open("x") as stream:
        for name, sha in sorted(hashes.items()):
            stream.write(f"{sha}  {name}\n")
    manifest = {
        "run": str(DEST), "binary_sha256": SOURCE_SHA,
        "accepted_checkpoint_sha256": CHECKPOINT_SHA,
        "accepted_time": .045616533734626756, "target_time": .05,
        "input_sha256": hashes,
        "scope": "Check selective preclip preserves historical corner-overflow recovery",
    }
    (DEST / "staging.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

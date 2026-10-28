"""Prestage a conditional t=6.5→10 native jitter control."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


DOCS = Path(__file__).resolve().parent
DEST = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/jitter005_to10_a10_v1")
PARENT_JOB = 437327
RUNNER = "kepler_q2_bounded_jitter005_to10_a10_20261027.slurm"
FILES = (
    RUNNER,
    "prepare_kepler_bounded_followon_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)
    DEST.mkdir()
    for name in FILES:
        shutil.copy2(DOCS / name, DEST / name)
    manifest = {
        "destination": str(DEST),
        "parent_job": PARENT_JOB,
        "dependency": f"afterok:{PARENT_JOB}",
        "runner": RUNNER,
        "prestaged_sha256": {name: digest(DEST / name) for name in FILES},
        "physics_state": "Copied only after the t=6.5 parent native audit passes",
    }
    (DEST / "prestage_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

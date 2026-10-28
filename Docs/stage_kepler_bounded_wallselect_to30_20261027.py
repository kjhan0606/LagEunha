"""Prestage a conditional t=20→30 fast CPU/selective CUDA continuation."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


DOCS = Path(__file__).resolve().parent
DEST = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to30_a10_v1")
PARENT_JOB = 437343
RUNNER = "kepler_q2_bounded_wallselect_to30_a10_20261027.slurm"
PREPARER = "prepare_kepler_bounded_wallselect_followon_20261027.py"
SOURCE = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/"
              "wallselect_pair_t58_a10_v1/candidate/eunha2.exe")
SOURCE_SHA = "87d60c1cb059e1710dd2749c84af5c98ba17fe077e6fe080b9c8b42d5fd40968"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    if digest(SOURCE) != SOURCE_SHA:
        raise ValueError("Selective CUDA binary changed")
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)
    DEST.mkdir()
    for name in (RUNNER, PREPARER):
        shutil.copy2(DOCS / name, DEST / name)
    manifest = {
        "destination": str(DEST), "parent_job": PARENT_JOB,
        "dependency": f"afterok:{PARENT_JOB}",
        "runner": RUNNER,
        "candidate_binary": str(SOURCE), "candidate_binary_sha256": SOURCE_SHA,
        "prestaged_sha256": {name: digest(DEST / name) for name in (RUNNER, PREPARER)},
        "physics_state": "Copied only after parent t=20 native physics audit passes",
    }
    (DEST / "prestage_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

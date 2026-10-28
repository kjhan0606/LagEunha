"""Stage a short replay from the accepted t=0.0456165 bounded checkpoint."""

import hashlib
import json
from pathlib import Path
import shutil


DOCS = Path("/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs")
SOURCE = Path("/tmp/Kepler_sparse_fixedpoly_20261027/eunha2.exe")
SOURCE_SHA = "8a038be1ab7b6b13bc55b3560fa2a88c607e195fed0b8d6a006b6a5f3b745872"
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
CHECKPOINT = ROOT / "to0825_v2/khdiag.step000024.start_geometry.dat"
CHECKPOINT_SHA = "08a7500e5d72c965c771ad51fd90d87a3cd7872a439edfb146df7a566d18c76a"
PARAMS = ROOT / "t001_v2/params.dat"
DEST = ROOT / "capacity_replay005_coarsefinal_v2"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    if digest(SOURCE) != SOURCE_SHA or digest(CHECKPOINT) != CHECKPOINT_SHA:
        raise ValueError("Binary or accepted checkpoint changed")
    files = {
        "eunha2.exe": SOURCE,
        "params.dat": PARAMS,
        "khout.000001.dat": CHECKPOINT,
        "job.slurm": DOCS / "kepler_q2_bounded_capacity_replay_to005_a10_20261027.slurm",
        "audit_kepler_bounded_short_replay_20261027.py":
            DOCS / "audit_kepler_bounded_short_replay_20261027.py",
        "audit_kepler_sparse_area_t0825_20261027.py":
            DOCS / "audit_kepler_sparse_area_t0825_20261027.py",
        "audit_kepler_sparse_state_consistency_20261027.py":
            DOCS / "audit_kepler_sparse_state_consistency_20261027.py",
        "prepare_kepler_sparse_area_restart_20261027.py":
            DOCS / "prepare_kepler_sparse_area_restart_20261027.py",
    }
    original_hashes = {name: digest(path) for name, path in files.items()}
    DEST.mkdir(exist_ok=False)
    for name, path in files.items():
        shutil.copy2(path, DEST/name)
        if digest(DEST/name) != original_hashes[name]:
            raise IOError(f"Staged {name} changed")
    (DEST/"inputs.sha256").write_text("".join(
        f"{value}  {name}\n" for name, value in sorted(original_hashes.items())))
    receipt = {"run": str(DEST), "binary_sha256": SOURCE_SHA,
               "accepted_checkpoint_sha256": CHECKPOINT_SHA,
               "accepted_time": 0.045616533734626756,
               "target_time": 0.05, "input_sha256": original_hashes}
    (DEST/"staging.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()

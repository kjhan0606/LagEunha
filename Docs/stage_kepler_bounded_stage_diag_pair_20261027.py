"""Stage paired RK4 stage-diagnostic timer runs from one accepted state."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import snapshot_time


DOCS = Path(__file__).resolve().parent
BASE = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/to0825_near32_v1")
PROFILE = Path("/tmp/Kepler_sparse_fixedpoly_profile_20261027/eunha2.exe")
CONTROL = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/rk_accel_cache_pair_t085_a10_v1/profile/eunha2.exe")
CONTROL_SHA = "58fa73d1b7f2ad36b478710b8ec4d0f2551300fe0e06ff530edb5c9ce664171b"
PROFILE_SHA = "30d81c891dfaa394d67f37ea6069eefd46ff93eddd4f46286e2056d7c5cd3626"
CHECKPOINT_SHA = "134e4e76a475c6923f1106ffb944fd0d51a0df82e7d4c86ea2814057ce9a12f6"
HELPERS = (
    "audit_kepler_bounded_short_replay_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)
RUNNER = "kepler_q2_bounded_stage_diag_pair_t085_a10_20261027.slurm"
PAIR_AUDIT = "audit_kepler_bounded_rk_gap_pair_20261027.py"


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
    if digest(CONTROL) != CONTROL_SHA:
        raise ValueError("Control binary changed")
    if digest(PROFILE) != PROFILE_SHA:
        raise ValueError("Profile binary changed")
    if (digest(source_snapshot) != CHECKPOINT_SHA or
            abs(snapshot_time(source_snapshot, 55351) - .825) > 1e-10):
        raise ValueError("Accepted checkpoint changed")
    receipt = json.loads((BASE / "audit_t0825.json").read_text())
    if receipt["endpoint_sha256"] != CHECKPOINT_SHA or not all(receipt["checks"].values()):
        raise ValueError("Predecessor physics gate failed")
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)
    dest.mkdir(parents=True)
    shutil.copy2(DOCS / RUNNER, dest / RUNNER)
    shutil.copy2(DOCS / PAIR_AUDIT, dest / PAIR_AUDIT)
    hashes = {}
    for variant, executable in (("control", CONTROL),
                                ("profile", PROFILE)):
        directory = dest / variant
        directory.mkdir()
        sources = {
            "eunha2.exe": executable,
            "params.dat": BASE / "params.dat",
            "khout.000001.dat": source_snapshot,
            **{name: DOCS / name for name in HELPERS},
        }
        hashes[variant] = {}
        for name, source in sources.items():
            shutil.copy2(source, directory / name)
            hashes[variant][name] = digest(directory / name)
        with (directory / "inputs.sha256").open("x") as stream:
            for name, sha in sorted(hashes[variant].items()):
                stream.write(f"{sha}  {name}\n")
    manifest = {
        "destination": str(dest), "source": str(BASE),
        "checkpoint_time": .825, "target_time": .85,
        "checkpoint_sha256": CHECKPOINT_SHA,
        "binaries": {"control": CONTROL_SHA, "profile": PROFILE_SHA},
        "runner_sha256": digest(dest / RUNNER),
        "pair_auditor_sha256": digest(dest / PAIR_AUDIT),
        "inputs": hashes,
        "scope": "Two sequential short replays on one A10; only RK4 stage-diagnostic timing differs",
    }
    (dest / "stage_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

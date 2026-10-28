"""Capture one native raw-face call for both selective-preclip variants."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PAIR = ROOT / "wallselect_pair_t085_a10_v1"
BASE = ROOT / "to0825_near32_v1"
DEST = ROOT / "wallselect_face_call1_a10_v1"
CHECKPOINT_SHA = "134e4e76a475c6923f1106ffb944fd0d51a0df82e7d4c86ea2814057ce9a12f6"
BINARIES = {
    "control": "56d96440005d82b424f00041f0631540ab22628932b238df285a54216ba4bc7b",
    "candidate": "87d60c1cb059e1710dd2749c84af5c98ba17fe077e6fe080b9c8b42d5fd40968",
}
RUNNER = "kepler_q2_bounded_wallselect_face_capture_a10_20261027.slurm"
HELPERS = (
    "audit_kepler_bounded_short_replay_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
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
    checkpoint = BASE / "khout.000403.dat"
    if digest(checkpoint) != CHECKPOINT_SHA or abs(snapshot_time(checkpoint, 55351)-.825) > 1e-10:
        raise ValueError("Wrong accepted checkpoint")
    receipt = json.loads((PAIR / "pair_audit.json").read_text())
    if not all(row["bitwise_equal"] for row in receipt["physical_fields"].values()):
        raise ValueError("Predecessor pair was not physically identical")
    for variant, sha in BINARIES.items():
        if digest(PAIR / variant / "eunha2.exe") != sha:
            raise ValueError(f"Changed {variant} binary")
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)
    DEST.mkdir()
    shutil.copy2(DOCS / RUNNER, DEST / RUNNER)
    staged = {}
    for variant in BINARIES:
        out = DEST / variant
        out.mkdir()
        files = {
            "eunha2.exe": PAIR / variant / "eunha2.exe",
            "params.dat": BASE / "params.dat",
            "khout.000001.dat": checkpoint,
            **{name: DOCS / name for name in HELPERS},
        }
        staged[variant] = {}
        for name, source in files.items():
            shutil.copy2(source, out / name)
            staged[variant][name] = digest(out / name)
        with (out / "inputs.sha256").open("x") as stream:
            for name, sha in sorted(staged[variant].items()):
                stream.write(f"{sha}  {name}\n")
    manifest = {
        "run": str(DEST), "checkpoint_sha256": CHECKPOINT_SHA,
        "start_time": .825, "stop_time": .826, "binary_sha256": BINARIES,
        "runner_sha256": digest(DEST / RUNNER), "inputs": staged,
        "scope": "One-force-call raw-face comparison; no long-time physics claim",
    }
    (DEST / "staging.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

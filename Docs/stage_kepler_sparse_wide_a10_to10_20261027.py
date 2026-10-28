"""Stage the accepted wide sparse Kepler t=10 continuation for one A10."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


FILES = (
    "kepler_q2_sparse_area_wide_to10_a10_20261027.slurm",
    "prepare_kepler_sparse_area_restart_20261027.py",
    "prepare_kepler_sparse_area_from0825_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_area_t10_20261027.py",
)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("destination", type=Path)
    args = ap.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(f"Refusing existing stage: {dest}")
    source = Path(__file__).resolve().parent
    dest.mkdir(parents=True)
    hashes = {}
    for name in FILES:
        shutil.copy2(source / name, dest / name)
        hashes[name] = digest(dest / name)
    with (dest / "scripts.sha256").open("x") as stream:
        for name, value in hashes.items():
            stream.write(f"{value}  {name}\n")
    subprocess.run(["bash", "-n", FILES[0]], cwd=dest, check=True)
    subprocess.run([sys.executable, "-c",
                    "import prepare_kepler_sparse_area_from0825_20261027; "
                    "import audit_kepler_sparse_area_t10_20261027"],
                   cwd=dest, check=True)
    manifest = {
        "destination": str(dest),
        "files_sha256": hashes,
        "validation": "bash syntax and Python imports passed",
        "scope": "Scripts only; native accepted t=.825 checkpoint and binary are copied by the runner after gates pass",
    }
    with (dest / "stage_manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

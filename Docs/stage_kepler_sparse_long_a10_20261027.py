"""Stage an immutable sparse Kepler long-run script bundle before sbatch."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


FILES = (
    "kepler_sparse_area_long_a10_20261027.slurm",
    "prepare_kepler_sparse_long_restart_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
    "audit_kepler_sparse_area_t10_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "audit_kepler_sparse_native_ballistic_20261027.py",
    "audit_kepler_seam_provenance_20261027.py",
    "audit_kepler_native_radial_20261026.py",
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("destination", type=Path)
    args = ap.parse_args()
    destination = args.destination.resolve()
    if destination.exists():
        raise ValueError("Refusing to stage over an existing directory")
    source = Path(__file__).resolve().parent
    destination.mkdir(parents=True)
    hashes = {}
    for name in FILES:
        shutil.copy2(source / name, destination / name)
        hashes[name] = sha256(destination / name)
    with (destination / "scripts.sha256").open("x") as stream:
        for name in FILES:
            stream.write(f"{hashes[name]}  {name}\n")
    subprocess.run([sys.executable, "-c",
                    "import prepare_kepler_sparse_long_restart_20261027; "
                    "import audit_kepler_sparse_area_t10_20261027; "
                    "import audit_kepler_sparse_state_consistency_20261027; "
                    "import audit_kepler_sparse_native_ballistic_20261027; "
                    "import audit_kepler_seam_provenance_20261027"],
                   cwd=destination, check=True)
    manifest = {"destination": str(destination), "files_sha256": hashes,
                "import_check": "passed",
                "scope": "A10 scripts only; native restart is copied after predecessor physics and state gates pass"}
    with (destination / "stage_manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

"""Stage one immutable audited H200 q2 segment without copying a checkpoint."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


FILES = (
    "kepler_polar_q2_halo28_segment_20261027.slurm",
    "prepare_kepler_polar_long_restart_20261027.py",
    "audit_kepler_polar_long_gate_20261027.py",
    "audit_kepler_native_radial_20261026.py",
    "audit_kepler_seam_checkpoint_20261026.py",
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
    source = Path(__file__).resolve().parent
    destination = args.destination.resolve()
    if destination.exists():
        raise ValueError(f"Refusing existing stage directory: {destination}")
    if not destination.parent.is_dir():
        raise ValueError("Stage parent must already exist")
    destination.mkdir()
    hashes = {}
    for name in FILES:
        shutil.copy2(source / name, destination / name)
        hashes[name] = digest(destination / name)
    with (destination / "scripts.sha256").open("x") as stream:
        for name in FILES:
            stream.write(f"{hashes[name]}  {name}\n")
    subprocess.run(["bash", "-n", FILES[0]], cwd=destination, check=True)
    subprocess.run(["sha256sum", "-c", "scripts.sha256"], cwd=destination,
                   stdout=subprocess.DEVNULL, check=True)
    subprocess.run([sys.executable, "-c",
                    "import prepare_kepler_polar_long_restart_20261027; "
                    "import audit_kepler_polar_long_gate_20261027; "
                    "import audit_kepler_native_radial_20261026"],
                   cwd=destination, check=True)
    manifest = {"destination": str(destination), "files_sha256": hashes,
                "bash_syntax": "passed", "sha256_check": "passed",
                "import_check": "passed",
                "scope": "Scripts only; restart copied by runner after parent physical gate"}
    with (destination / "stage_manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

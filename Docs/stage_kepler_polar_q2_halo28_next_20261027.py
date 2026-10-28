"""Stage the next guarded two-unit q2 W0 continuation segment."""

import argparse
import hashlib
import json
from decimal import Decimal
from pathlib import Path
import re
import shutil
import subprocess
import sys


EXPECTED_PARAMS = "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"
EXPECTED_BINARY = "7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284"
FILES = (
    "job.slurm",
    "prepare_kepler_polar_long_restart_20261027.py",
    "audit_kepler_polar_long_gate_20261027.py",
    "audit_kepler_native_radial_20261026.py",
    "audit_kepler_seam_checkpoint_20261026.py",
)
CASE_RE = re.compile(
    r"^predecessor_time=([^;]+); target_time=([^;]+); "
    r"predecessor_receipt=([^;]+); tag=([^\s]+)$",
    re.MULTILINE,
)


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_hashes(path):
    result = {}
    for line in path.read_text().splitlines():
        sha256, name = line.split(maxsplit=1)
        result[name.strip()] = sha256
    return result


def label(value):
    number = Decimal(str(value))
    if number == number.to_integral_value():
        return str(int(number))
    return format(number.normalize(), "f").replace(".", "p")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("parent", type=Path)
    parser.add_argument("dependency_job", type=int)
    parser.add_argument("stop", type=Decimal)
    args = parser.parse_args()

    parent = args.parent.resolve()
    start = args.stop - Decimal("2")
    if args.stop <= 20 or args.stop != args.stop.to_integral_value():
        raise ValueError("This stager expects an integer endpoint after t=20")
    if not parent.is_dir():
        raise FileNotFoundError(parent)

    parent_hashes = read_hashes(parent / "scripts.sha256")
    for name, expected in parent_hashes.items():
        if digest(parent / name) != expected:
            raise ValueError(f"Parent stage input changed: {name}")

    parent_manifest = json.loads((parent / "stage_manifest.json").read_text())
    if parent_manifest.get("destination") != str(parent):
        raise ValueError("Unexpected predecessor manifest destination")
    if parent_manifest.get("files_sha256") != parent_hashes:
        raise ValueError("Predecessor manifest and scripts.sha256 disagree")
    if parent_manifest.get("expected_params_sha256") != EXPECTED_PARAMS:
        raise ValueError("Unexpected predecessor parameter hash")
    if parent_manifest.get("expected_binary_sha256") != EXPECTED_BINARY:
        raise ValueError("Unexpected predecessor executable hash")

    job = (parent / "job.slurm").read_text()
    match = CASE_RE.search(job)
    if not match:
        raise ValueError("Could not parse predecessor segment contract")
    prev_start, prev_stop, prev_receipt, prev_tag = match.groups()
    if Decimal(prev_stop) != start:
        raise ValueError("Requested continuation does not start at parent endpoint")
    start_label = label(start)
    stop_label = label(args.stop)
    parent_start_label = label(Decimal(prev_start))
    if (prev_receipt != f"audit_t{parent_start_label}.json"
            or prev_tag != f"t{start_label}"):
        raise ValueError("Predecessor input receipt or endpoint tag is inconsistent")
    if job.count(f"#SBATCH --job-name=LagQ2H28To{start_label}") != 1:
        raise ValueError("Unexpected predecessor job name")

    old_case = match.group(0)
    new_case = (
        f"predecessor_time={start_label}; target_time={stop_label}; "
        f"predecessor_receipt=audit_t{start_label}.json; tag=t{stop_label}"
    )
    job = job.replace(old_case, new_case, 1)
    job = job.replace(f"#SBATCH --job-name=LagQ2H28To{start_label}",
                      f"#SBATCH --job-name=LagQ2H28To{stop_label}", 1)

    run = parent / f"to{stop_label}_h28_guarded_20261027"
    temp = parent / f"to{stop_label}_h28_guarded_20261027.stage_tmp"
    if run.exists() or temp.exists():
        raise FileExistsError("Refusing an existing continuation or staging directory")

    temp.mkdir()
    try:
        hashes = {}
        for name in FILES:
            destination = temp / name
            if name == "job.slurm":
                destination.write_text(job)
            else:
                shutil.copy2(parent / name, destination)
                if digest(destination) != parent_hashes.get(name):
                    raise ValueError(f"Changed validated helper: {name}")
            hashes[name] = digest(destination)

        subprocess.run(["bash", "-n", "job.slurm"], cwd=temp, check=True)
        with (temp / "scripts.sha256").open("x") as stream:
            for name in FILES:
                stream.write(f"{hashes[name]}  {name}\n")
        subprocess.run(["sha256sum", "-c", "scripts.sha256"], cwd=temp,
                       stdout=subprocess.DEVNULL, check=True)
        subprocess.run([sys.executable, "-c",
                        "import prepare_kepler_polar_long_restart_20261027; "
                        "import audit_kepler_polar_long_gate_20261027; "
                        "import audit_kepler_native_radial_20261026"],
                       cwd=temp, check=True)

        manifest = {
            "destination": str(run),
            "parent": str(parent),
            "predecessor_receipt": f"audit_t{start_label}.json",
            "start": float(start),
            "stop": float(args.stop),
            "dependency_job": args.dependency_job,
            "expected_params_sha256": EXPECTED_PARAMS,
            "expected_binary_sha256": EXPECTED_BINARY,
            "files_sha256": hashes,
            "bash_syntax": "passed",
            "sha256_check": "passed",
            "import_check": "passed",
            "scope": "Full-domain W0/Voronoi, fixed-mass q2 continuation; restart is gated by the predecessor physics receipt.",
        }
        with (temp / "stage_manifest.json").open("x") as stream:
            json.dump(manifest, stream, indent=2)
            stream.write("\n")
        temp.rename(run)
    except Exception:
        shutil.rmtree(temp, ignore_errors=True)
        raise

    print(json.dumps({"run": str(run), "start": float(start),
                      "stop": float(args.stop),
                      "dependency_job": args.dependency_job,
                      "files_sha256": hashes,
                      "validation": "passed"}, indent=2))


if __name__ == "__main__":
    main()

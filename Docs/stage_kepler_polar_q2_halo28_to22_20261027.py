"""Stage the next guarded q2 W0 interval, t=20 to t=22."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


PARENT = Path(
    "/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/"
    "halo12_after_failure_to0825/to2p221_h28/"
    "to4p5_h28_guarded_20261027/to6p5_h28_guarded_20261027/"
    "to8p5_h28_guarded_20261027/to10p0_h28_guarded_20261027/"
    "to12_h28_guarded_20261027/to14_h28_guarded_20261027/"
    "to16_h28_guarded_20261027/to18_h28_guarded_20261027/"
    "to20_h28_guarded_20261027"
)
RUN = PARENT / "to22_h28_guarded_20261027"
TEMP = PARENT / "to22_h28_guarded_20261027.stage_tmp"
DEPENDENCY_JOB = 437874
EXPECTED_PARAMS = "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"
EXPECTED_BINARY = "7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284"
FILES = (
    "job.slurm",
    "prepare_kepler_polar_long_restart_20261027.py",
    "audit_kepler_polar_long_gate_20261027.py",
    "audit_kepler_native_radial_20261026.py",
    "audit_kepler_seam_checkpoint_20261026.py",
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


def main():
    if RUN.exists() or TEMP.exists():
        raise FileExistsError("Refusing an existing t=22 stage or staging directory")
    if not PARENT.is_dir():
        raise FileNotFoundError(PARENT)

    parent_hashes = read_hashes(PARENT / "scripts.sha256")
    for name, expected in parent_hashes.items():
        if digest(PARENT / name) != expected:
            raise ValueError(f"Parent stage input changed: {name}")

    parent_manifest = json.loads((PARENT / "stage_manifest.json").read_text())
    if parent_manifest.get("destination") != str(PARENT):
        raise ValueError("Unexpected t=20 parent manifest")
    if parent_manifest.get("files_sha256") != parent_hashes:
        raise ValueError("t=20 manifest and scripts.sha256 disagree")

    job = (PARENT / "job.slurm").read_text()
    old_case = (
        "predecessor_time=18; target_time=20; "
        "predecessor_receipt=audit_t18.json; tag=t20"
    )
    new_case = (
        "predecessor_time=20; target_time=22.0; "
        "predecessor_receipt=audit_t20.json; tag=t22"
    )
    if job.count(old_case) != 1 or job.count("#SBATCH --job-name=LagQ2H28To20") != 1:
        raise ValueError("Unexpected t=20 runner template")
    job = job.replace(old_case, new_case, 1)
    job = job.replace("#SBATCH --job-name=LagQ2H28To20",
                      "#SBATCH --job-name=LagQ2H28To22", 1)

    TEMP.mkdir()
    try:
        hashes = {}
        for name in FILES:
            destination = TEMP / name
            if name == "job.slurm":
                destination.write_text(job)
            else:
                shutil.copy2(PARENT / name, destination)
                if digest(destination) != parent_hashes.get(name):
                    raise ValueError(f"Changed validated helper: {name}")
            hashes[name] = digest(destination)

        subprocess.run(["bash", "-n", "job.slurm"], cwd=TEMP, check=True)
        with (TEMP / "scripts.sha256").open("x") as stream:
            for name in FILES:
                stream.write(f"{hashes[name]}  {name}\n")
        subprocess.run(["sha256sum", "-c", "scripts.sha256"], cwd=TEMP,
                       stdout=subprocess.DEVNULL, check=True)
        subprocess.run([sys.executable, "-c",
                        "import prepare_kepler_polar_long_restart_20261027; "
                        "import audit_kepler_polar_long_gate_20261027; "
                        "import audit_kepler_native_radial_20261026"],
                       cwd=TEMP, check=True)

        manifest = {
            "destination": str(RUN),
            "parent": str(PARENT),
            "predecessor_receipt": "audit_t20.json",
            "start": 20.0,
            "stop": 22.0,
            "dependency_job": DEPENDENCY_JOB,
            "expected_params_sha256": EXPECTED_PARAMS,
            "expected_binary_sha256": EXPECTED_BINARY,
            "files_sha256": hashes,
            "bash_syntax": "passed",
            "sha256_check": "passed",
            "import_check": "passed",
            "scope": "Full-domain W0/Voronoi, fixed-mass q2 continuation; restart is gated by the t=20 physical receipt.",
        }
        with (TEMP / "stage_manifest.json").open("x") as stream:
            json.dump(manifest, stream, indent=2)
            stream.write("\n")
        TEMP.rename(RUN)
    except Exception:
        shutil.rmtree(TEMP, ignore_errors=True)
        raise

    print(json.dumps({"run": str(RUN), "start": 20.0, "stop": 22.0,
                      "dependency_job": DEPENDENCY_JOB,
                      "files_sha256": hashes,
                      "validation": "passed"}, indent=2))


if __name__ == "__main__":
    main()

"""Stage the guarded full-domain W0, halo-28 continuation from t=6.5 to 10."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path("/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825")
Q2_ROOT = Path("/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2")
PARENT = ROOT / "to2p221_h28/to4p5_h28_guarded_20261027/to6p5_h28_guarded_20261027"
RUN = PARENT / "to10_h28_guarded_20261027"
DOCS = Path(__file__).resolve().parent
DEPENDENCY_JOB = 437251
EXPECTED_PARAMS = "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"
EXPECTED_BINARY = "7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284"
INITIAL = Q2_ROOT / "khout.000000.dat"
FILES = (
    "job.slurm",
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


def read_hashes(path):
    result = {}
    for line in path.read_text().splitlines():
        sha256, name = line.split(maxsplit=1)
        result[name.strip()] = sha256
    return result


def main():
    if RUN.exists():
        raise FileExistsError(RUN)
    if not PARENT.is_dir() or not INITIAL.is_file():
        raise FileNotFoundError("Expected t=6.5 stage parent or polar initial snapshot")
    parent_scripts = read_hashes(PARENT / "scripts.sha256")
    for name, expected in parent_scripts.items():
        if digest(PARENT / name) != expected:
            raise ValueError(f"Parent stage script changed: {name}")
    if parent_scripts.get("prepare_kepler_polar_long_restart_20261027.py") is None:
        raise ValueError("Missing validated restart helper in t=6.5 stage")

    job = (PARENT / "kepler_polar_q2_halo28_segment_20261027.slurm").read_text()
    if job != (DOCS / "kepler_polar_q2_halo28_segment_20261027.slurm").read_text():
        raise ValueError("Validated H28 segment template differs from Docs")
    if job.count("#SBATCH --job-name=LagQ2H28Segment") != 1:
        raise ValueError("Unexpected H28 job-name template")
    start = job.index('if [ "$#" -ne 1 ]; then')
    end = job.index("\nsha256sum -c scripts.sha256", start)
    fixed_segment = (
        "predecessor_time=6.5; target_time=10.0; "
        "predecessor_receipt=audit_t6p5.json; tag=t10\n"
    )
    job = (job[:start] + fixed_segment + job[end + 1:])
    job = job.replace("#SBATCH --job-name=LagQ2H28Segment",
                      "#SBATCH --job-name=LagQ2H28To10", 1)

    RUN.mkdir()
    hashes = {}
    for name in FILES:
        if name == "job.slurm":
            destination = RUN / name
            destination.write_text(job)
            hashes[name] = digest(destination)
        else:
            source = PARENT / name
            shutil.copy2(source, RUN / name)
            hashes[name] = digest(RUN / name)
            if hashes[name] != parent_scripts.get(name):
                raise ValueError(f"Changed validated helper: {name}")

    subprocess.run(["bash", "-n", "job.slurm"], cwd=RUN, check=True)
    with (RUN / "scripts.sha256").open("x") as stream:
        for name in FILES:
            stream.write(f"{hashes[name]}  {name}\n")
    subprocess.run(["sha256sum", "-c", "scripts.sha256"], cwd=RUN,
                   stdout=subprocess.DEVNULL, check=True)
    with (RUN / "stage_manifest.json").open("x") as stream:
        json.dump({"destination": str(RUN), "parent": str(PARENT),
                   "predecessor_receipt": "audit_t6p5.json",
                   "start": 6.5, "stop": 10.0,
                   "dependency_job": DEPENDENCY_JOB,
                   "initial_sha256": digest(INITIAL),
                   "params_sha256": EXPECTED_PARAMS,
                   "expected_binary_sha256": EXPECTED_BINARY,
                   "files_sha256": hashes,
                   "bash_syntax": "passed", "sha256_check": "passed",
                   "scope": "Full-domain W0/Voronoi, halo-28 segment at unchanged disk-annulus resolution; restart follows the t=6.5 physical gate."},
                  stream, indent=2)
        stream.write("\n")
    print(json.dumps({"run": str(RUN), "start": 6.5, "stop": 10.0,
                      "dependency_job": DEPENDENCY_JOB,
                      "sha256_check": "passed"}, indent=2))


if __name__ == "__main__":
    main()

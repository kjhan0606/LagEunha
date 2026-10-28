"""Stage a guarded one-rank/two-rank v11 probe from the accepted t=40 state."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DOCS = Path(__file__).resolve().parent
PARENT = ROOT / "gpu_energy_neutral_thermalhoist_v11_fresh_to40_h200_v2"
RUN = ROOT / "gpu_energy_neutral_v11_mpi_rank_pair_t40_h200_v2"
EXPECTED_BINARY = "adeca5d738bce85a740b966fbbef04e3aaaba898ac7978eb5db797a620eaf5cd"
EXPECTED_INITIAL = "8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b"
DEPENDENCY_JOB = 437796


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


def copy_checked(source, destination, expected=None):
    if not source.is_file():
        raise FileNotFoundError(source)
    shutil.copy2(source, destination)
    actual = digest(destination)
    if actual != digest(source) or (expected is not None and actual != expected):
        raise ValueError(f"Copy/hash mismatch: {source} -> {destination}")
    return actual


def write_manifest(path, hashes):
    with path.open("x") as output:
        for name, value in sorted(hashes.items()):
            output.write(f"{value}  {name}\n")


def main():
    parent_manifest = json.loads((PARENT / "stage_manifest.json").read_text())
    if (parent_manifest.get("start") != 30.025
            or parent_manifest.get("stop") != 40.0
            or parent_manifest.get("binary_sha256") != EXPECTED_BINARY
            or parent_manifest.get("initial_sha256") != EXPECTED_INITIAL
            or parent_manifest.get("omp_threads") != 8):
        raise ValueError("Unexpected v11 t=30.025→40 parent lineage")

    parent_inputs = read_hashes(PARENT / "inputs.sha256")
    if (parent_inputs.get("eunha2.exe") != EXPECTED_BINARY
            or parent_inputs.get("initial.dat") != EXPECTED_INITIAL):
        raise ValueError("Parent input manifest changed")
    binary = PARENT / "eunha2.exe"
    if digest(binary) != EXPECTED_BINARY:
        raise ValueError("Parent v11 executable missing or changed")
    if RUN.exists():
        raise FileExistsError(RUN)

    job_source = DOCS / "kepler_weighted_v11_mpi_scaling_t40_h200_20261027.slurm"
    subprocess.run(["bash", "-n", str(job_source)], check=True)
    RUN.mkdir()

    root_sources = {
        "job.slurm": job_source,
        "audit_kepler_weighted_mpi_rank_pair_20261027.py":
            DOCS / "audit_kepler_weighted_mpi_rank_pair_20261027.py",
        "audit_kepler_sparse_area_t0825_20261027.py":
            PARENT / "audit_kepler_sparse_area_t0825_20261027.py",
        "prepare_kepler_sparse_area_restart_20261027.py":
            PARENT / "prepare_kepler_sparse_area_restart_20261027.py",
    }
    root_hashes = {}
    for name, source in root_sources.items():
        root_hashes[name] = copy_checked(source, RUN / name)
    write_manifest(RUN / "inputs.sha256", root_hashes)

    common_names = [name for name in parent_inputs
                    if name not in ("job.slurm", "eunha2.exe")]
    branch_hashes = {}
    for label in ("control", "candidate"):
        directory = RUN / label
        directory.mkdir()
        hashes = {}
        for name in common_names:
            hashes[name] = copy_checked(PARENT / name, directory / name,
                                        expected=parent_inputs[name])
        hashes["eunha2.exe"] = copy_checked(binary, directory / "eunha2.exe",
                                              expected=EXPECTED_BINARY)
        write_manifest(directory / "inputs.sha256", hashes)
        branch_hashes[label] = hashes

    with (RUN / "stage_manifest.json").open("x") as output:
        json.dump({"parent_run": str(PARENT), "start": 40.0,
                   "stop": 40.25, "dependency_job": DEPENDENCY_JOB,
                   "binary_sha256": EXPECTED_BINARY,
                   "initial_sha256": EXPECTED_INITIAL,
                   "ranks": {"control": 1, "candidate": 2},
                   "omp_threads_per_rank": 8,
                   "root_inputs": root_hashes,
                   "branch_inputs": branch_hashes}, output, indent=2)
        output.write("\n")
    print(RUN)


if __name__ == "__main__":
    main()

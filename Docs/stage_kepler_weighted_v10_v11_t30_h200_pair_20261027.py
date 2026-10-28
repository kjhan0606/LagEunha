"""Stage a matched v10/v11 native pair from the accepted weighted t=30 state."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DOCS = Path(__file__).resolve().parent
PARENT = ROOT / "gpu_energy_neutral_capacity64_v5_fresh_to30_omp8_h200_v1"
V10_BUILD = ROOT / "gpu_energy_neutral_wallselect_cache_capacity64_build_v10"
V11_BUILD = ROOT / "gpu_energy_neutral_wallselect_thermalhoist_capacity64_build_v11"
RUN = ROOT / "gpu_energy_neutral_thermalhoist_pair_t30_h200_v2"
EXPECTED_PARENT = "feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b"
EXPECTED_INITIAL = "8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b"


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
    parent_manifest = json.loads((PARENT / "stage_manifest.json").read_text())
    if (parent_manifest.get("start") != 20.0
            or parent_manifest.get("stop") != 30.0
            or parent_manifest.get("binary_sha256") != EXPECTED_PARENT
            or parent_manifest.get("initial_sha256") != EXPECTED_INITIAL
            or parent_manifest.get("omp_threads") != 8):
        raise ValueError("Unexpected accepted t=20→30 parent lineage")
    parent_inputs = read_hashes(PARENT / "inputs.sha256")
    if (parent_inputs.get("eunha2.exe") != EXPECTED_PARENT
            or parent_inputs.get("initial.dat") != EXPECTED_INITIAL):
        raise ValueError("Parent input manifest changed")

    v10_sha = (V10_BUILD / "binary.sha256").read_text().split()[0]
    v11_sha = (V11_BUILD / "binary.sha256").read_text().split()[0]
    v10 = V10_BUILD / "eunha2_wallselect_cache_capacity64_v10.exe"
    v11 = V11_BUILD / "eunha2_wallselect_cache_thermal_capacity64_v11.exe"
    if (v10_sha != "0e92efe71e2516b8a2f1b919ea6ba1c5c019532d1dd407331988b996dda46775"
            or digest(v10) != v10_sha):
        raise ValueError("v10 binary missing or changed")
    if (v11_sha != "adeca5d738bce85a740b966fbbef04e3aaaba898ac7978eb5db797a620eaf5cd"
            or digest(v11) != v11_sha):
        raise ValueError("v11 binary missing or changed")
    if RUN.exists():
        raise FileExistsError(RUN)

    job_source = DOCS / "kepler_weighted_v10_v11_t30_h200_pair_20261027.slurm"
    subprocess.run(["bash", "-n", str(job_source)], check=True)
    RUN.mkdir()
    root_sources = {
        "job.slurm": job_source,
        "audit_kepler_weighted_fast_pair_20261027.py":
            DOCS / "audit_kepler_weighted_fast_pair_20261027.py",
        "prepare_kepler_sparse_area_restart_20261027.py":
            PARENT / "prepare_kepler_sparse_area_restart_20261027.py",
    }
    root_hashes = {}
    for name, source in root_sources.items():
        destination = RUN / name
        shutil.copy2(source, destination)
        root_hashes[name] = digest(destination)
        if root_hashes[name] != digest(source):
            raise ValueError(f"Changed root input: {name}")
    with (RUN / "inputs.sha256").open("x") as output:
        for name, value in sorted(root_hashes.items()):
            output.write(f"{value}  {name}\n")

    common_names = [name for name in parent_inputs
                    if name not in ("job.slurm", "eunha2.exe")]
    branch_hashes = {}
    for label, binary in (("control", v10), ("candidate", v11)):
        directory = RUN / label
        directory.mkdir()
        sources = {name: PARENT / name for name in common_names}
        sources["eunha2.exe"] = binary
        hashes = {}
        for name, source in sources.items():
            if not source.is_file():
                raise FileNotFoundError(source)
            destination = directory / name
            shutil.copy2(source, destination)
            hashes[name] = digest(destination)
            expected = parent_inputs.get(name)
            if name != "eunha2.exe" and expected and hashes[name] != expected:
                raise ValueError(f"Parent input mismatch: {label}/{name}")
            if hashes[name] != digest(source):
                raise ValueError(f"Changed branch input: {label}/{name}")
        with (directory / "inputs.sha256").open("x") as output:
            for name, value in sorted(hashes.items()):
                output.write(f"{value}  {name}\n")
        branch_hashes[label] = hashes

    with (RUN / "stage_manifest.json").open("x") as output:
        json.dump({"parent_run": str(PARENT), "parent_time": 30.0,
                   "target_time": 30.025, "dependency_job": 437613,
                   "initial_sha256": EXPECTED_INITIAL,
                   "parent_binary_sha256": EXPECTED_PARENT,
                   "control_binary_sha256": v10_sha,
                   "candidate_binary_sha256": v11_sha,
                   "root_inputs": root_hashes,
                   "branch_inputs": branch_hashes}, output, indent=2)
        output.write("\n")
    print(RUN)


if __name__ == "__main__":
    main()

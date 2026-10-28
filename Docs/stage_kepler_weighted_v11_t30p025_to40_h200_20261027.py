"""Stage the guarded v11 continuation after the matched t=30 pair."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DOCS = Path(__file__).resolve().parent
PAIR = ROOT / "gpu_energy_neutral_thermalhoist_pair_t30_h200_v2"
BUILD = ROOT / "gpu_energy_neutral_wallselect_thermalhoist_capacity64_build_v11"
RUN = ROOT / "gpu_energy_neutral_thermalhoist_v11_fresh_to40_h200_v2"
EXPECTED_BINARY = "adeca5d738bce85a740b966fbbef04e3aaaba898ac7978eb5db797a620eaf5cd"
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
    binary = BUILD / "eunha2_wallselect_cache_thermal_capacity64_v11.exe"
    binary_sha = (BUILD / "binary.sha256").read_text().split()[0]
    if binary_sha != EXPECTED_BINARY or digest(binary) != EXPECTED_BINARY:
        raise ValueError("v11 binary missing or changed")
    pair_manifest = json.loads((PAIR / "stage_manifest.json").read_text())
    if (pair_manifest.get("parent_time") != 30.0
            or pair_manifest.get("target_time") != 30.025
            or pair_manifest.get("candidate_binary_sha256") != EXPECTED_BINARY
            or pair_manifest.get("initial_sha256") != EXPECTED_INITIAL):
        raise ValueError("Unexpected v10/v11 pair lineage")
    candidate_dir = PAIR / "candidate"
    parent_hashes = read_hashes(candidate_dir / "inputs.sha256")
    if (parent_hashes.get("eunha2.exe") != EXPECTED_BINARY
            or parent_hashes.get("initial.dat") != EXPECTED_INITIAL):
        raise ValueError("Pair candidate input manifest changed")
    if RUN.exists():
        raise FileExistsError(RUN)

    job_source = DOCS / "kepler_weighted_v11_t30p025_to40_h200_20261027.slurm"
    subprocess.run(["bash", "-n", str(job_source)], check=True)
    RUN.mkdir()
    sources = {name: candidate_dir / name for name in parent_hashes
               if name not in ("eunha2.exe", "khout.000001.dat")}
    sources["eunha2.exe"] = binary
    sources["job.slurm"] = job_source
    hashes = {}
    for name, source in sources.items():
        destination = RUN / name
        shutil.copy2(source, destination)
        hashes[name] = digest(destination)
        expected = parent_hashes.get(name)
        if name not in ("eunha2.exe", "job.slurm") and expected and hashes[name] != expected:
            raise ValueError(f"Pair input mismatch: {name}")
        if hashes[name] != digest(source):
            raise ValueError(f"Changed staged input: {name}")
    with (RUN / "inputs.sha256").open("x") as output:
        for name, value in sorted(hashes.items()):
            output.write(f"{value}  {name}\n")
    with (RUN / "stage_manifest.json").open("x") as output:
        json.dump({"parent_pair": str(PAIR), "start": 30.025, "stop": 40.0,
                   "dependency_job": 437794, "binary_sha256": EXPECTED_BINARY,
                   "initial_sha256": EXPECTED_INITIAL, "omp_threads": 8,
                   "input_sha256": hashes}, output, indent=2)
        output.write("\n")
    print(RUN)


if __name__ == "__main__":
    main()

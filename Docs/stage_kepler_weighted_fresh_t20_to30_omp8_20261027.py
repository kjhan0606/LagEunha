"""Stage a t=20→30 weighted restart using eight OpenMP threads."""
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PARENT = ROOT / "gpu_energy_neutral_capacity64_v5_fresh_to20_a10_v1"
RUN = ROOT / "gpu_energy_neutral_capacity64_v5_fresh_to30_omp8_a10_v1"
FILES = (
    "kepler_fresh_weighted_t20_to30_omp8_20261027.slurm",
    "audit_kepler_bounded_native_evolution_20261027.py",
    "audit_kepler_q2_weighted_equilibrium_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_native_ballistic_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "disk_active_mask.bin",
    "eunha2.exe",
    "initial.dat",
    "laguerre_nonlinear_volume_probe_t6_20261026.py",
    "params.dat",
    "prepare_kepler_sparse_area_restart_20261027.py",
)
EXPECTED_BINARY = "feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b"
EXPECTED_INITIAL = "8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b"


def main():
    if RUN.exists():
        raise FileExistsError(RUN)
    parent_hashes = {}
    for line in (PARENT / "inputs.sha256").read_text().splitlines():
        sha256, name = line.split(maxsplit=1)
        parent_hashes[name.strip()] = sha256
    for name, expected in (("eunha2.exe", EXPECTED_BINARY),
                           ("initial.dat", EXPECTED_INITIAL)):
        if parent_hashes.get(name) != expected:
            raise ValueError(f"Unexpected t=20 parent input: {name}")
    subprocess.run(["bash", "-n", str(DOCS / FILES[0])], check=True)
    RUN.mkdir()
    hashes = {}
    for name in FILES:
        source = DOCS / name if name in FILES[:1] else PARENT / name
        destination = RUN / name
        shutil.copy2(source, destination)
        hashes[name] = digest(destination)
        if hashes[name] != digest(source):
            raise ValueError(f"Changed staged file: {name}")
        if name in parent_hashes and hashes[name] != parent_hashes[name]:
            raise ValueError(f"Parent input hash differs: {name}")
    with (RUN / "inputs.sha256").open("x") as output:
        for name, sha256 in sorted(hashes.items()):
            output.write(f"{sha256}  {name}\n")
    with (RUN / "stage_manifest.json").open("x") as output:
        json.dump({"parent_run": str(PARENT), "start": 20.0, "stop": 30.0,
                   "binary_sha256": hashes["eunha2.exe"], "initial_sha256": hashes["initial.dat"],
                   "omp_threads": 8, "input_sha256": hashes}, output, indent=2)
        output.write("\n")
    print(RUN)


if __name__ == "__main__":
    main()

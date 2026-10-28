"""Stage a Jacobi64 native branch from audited weighted t=7.426 snapshots."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PARENT = ROOT / "gpu_energy_neutral_weighted_t10_a10_v1"
BUILD = ROOT / "gpu_energy_neutral_native_build_v3_jacobi_choice"
BINARY = BUILD / "eunha2_energy_neutral_v3.exe"
RUNNER = "kepler_q2_weighted_jacobi64_from7426_to10_a10_20261027.slurm"
HELPERS = (
    "audit_kepler_q2_weighted_t10_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "laguerre_nonlinear_volume_probe_t6_20261026.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)
PARAMS_SHA = "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"


def copy_inputs(destination, sources):
    destination.mkdir()
    hashes = {}
    for name, source in sources.items():
        shutil.copy2(source, destination / name)
        hashes[name] = digest(destination / name)
    with (destination / "inputs.sha256").open("x") as output:
        for name, value in sorted(hashes.items()):
            output.write(f"{value}  {name}\n")
    return hashes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("destination", type=Path)
    args = ap.parse_args()
    destination = args.destination.resolve()
    if destination.exists():
        raise FileExistsError(destination)
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)
    build_receipt = (BUILD / "binary.sha256").read_text().split()[0]
    if digest(BINARY) != build_receipt:
        raise ValueError("Changed Jacobi-choice binary")
    monitor = json.loads((DOCS / "KEPLER_Q2_WEIGHTED_STEP000478_MONITOR_2026-10-27.json").read_text())
    jacobi = json.loads((DOCS / "KEPLER_Q2_WEIGHTED_JACOBI_RESIDUAL_T7426_2026-10-27.json").read_text())
    accepted = {}
    for label in ("aligned", "jitter"):
        source = PARENT / label
        checkpoint = source / "khout.000478.dat"
        initial = source / "initial.dat"
        t = snapshot_time(checkpoint, 55351)
        if (abs(t - 7.4262993674) > 1e-8
                or digest(checkpoint) != monitor["cases"][label]["weighted_endpoint_sha256"]
                or digest(checkpoint) != jacobi["cases"][label]["checkpoint_sha256"]
                or digest(initial) != monitor["cases"][label]["initial_sha256"]
                or digest(source / "params.dat") != PARAMS_SHA
                or jacobi["cases"][label]["independent_area_error_max"] >= 1e-6
                or jacobi["cases"][label]["projected_jacobi"]["64"]["jacobi_relative_residual_l2"] >= .1):
            raise ValueError(f"{label}: unaccepted or changed t=7.426 checkpoint")
        accepted[label] = (source, initial, checkpoint, t)
    destination.mkdir(parents=True)
    result = {
        "run": str(destination), "binary_sha256": build_receipt,
        "start_time_by_case": {}, "target_time": 10.,
        "jacobi_iterations": 64, "cases": {},
        "scope": "Branches from the same audited weighted t=7.426 checkpoints; initial t=0 field retained for endpoint metrics.",
    }
    for label, (source, initial, checkpoint, t) in accepted.items():
        case = destination / label
        hashes = copy_inputs(case, {
            "eunha2.exe": BINARY,
            "params.dat": source / "params.dat",
            "khout.000001.dat": checkpoint,
            "initial.dat": initial,
            "disk_active_mask.bin": source / "disk_active_mask.bin",
            "job.slurm": DOCS / RUNNER,
            **{name: DOCS / name for name in HELPERS},
        })
        result["start_time_by_case"][label] = t
        result["cases"][label] = hashes
    with (destination / "stage_manifest.json").open("x") as output:
        json.dump(result, output, indent=2)
        output.write("\n")
    print(json.dumps({"run": str(destination), "binary_sha256": build_receipt,
                      "start_time_by_case": result["start_time_by_case"]}, indent=2))


if __name__ == "__main__":
    main()

"""Stage independent aligned/jitter weighted t=6.501→10 native continuations."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PARENT = ROOT / "gpu_energy_neutral_weighted_t6501_guard_a10_v2"
BUILD = ROOT / "gpu_energy_neutral_guard_build_v2"
BINARY = BUILD / "eunha2_energy_neutral_guard.exe"
RUNNER = "kepler_q2_weighted_to10_a10_20261027.slurm"
HELPERS = (
    "audit_kepler_q2_weighted_long_20261027.py",
    "audit_kepler_bounded_native_evolution_20261027.py",
    "audit_kepler_sparse_native_ballistic_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "laguerre_nonlinear_volume_probe_t6_20261026.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)
INITIAL_RECEIPTS = {
    "aligned": "KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json",
    "jitter": "KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json",
}


def copy_inputs(destination, sources):
    destination.mkdir(exist_ok=True)
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
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)
    binary_sha = (BUILD / "binary.sha256").read_text().split()[0]
    if digest(BINARY) != binary_sha:
        raise ValueError("Native weighted binary changed")
    parent_manifest = json.loads((PARENT / "stage_manifest.json").read_text())
    mask = PARENT / "disk_active_mask.bin"
    if (mask.stat().st_size != 55351
            or digest(mask) != parent_manifest["mask_sha256"]
            or binary_sha != parent_manifest["binary_sha256"]):
        raise ValueError("Weighted mask or binary lineage changed")
    accepted = {}
    for label, receipt_name in INITIAL_RECEIPTS.items():
        one = json.loads((PARENT / label / "audit_t6501.json").read_text())
        initial = json.loads((DOCS / receipt_name).read_text())
        parent = Path(one["endpoint"])
        reference = Path(initial["initial"])
        if (not all(one["checks"].values()) or not all(initial["checks"].values())
                or parent.parent != PARENT / label
                or digest(parent) != one["endpoint_sha256"]
                or one["restart_sha256"] != initial["endpoint_sha256"]
                or digest(reference) != initial["initial_sha256"]
                or abs(snapshot_time(parent, 55351)-6.501)>1e-10
                or abs(snapshot_time(reference, 55351))>1e-12):
            raise ValueError(f"{label}: rejected or changed parent lineage")
        accepted[label] = (one, initial, parent, reference)
    dest.mkdir(parents=True)
    root_inputs = copy_inputs(dest, {
        "disk_active_mask.bin": mask,
        "job.slurm": DOCS / RUNNER,
    })
    cases = {}
    for label, (one, initial, parent, reference) in accepted.items():
        cases[label] = copy_inputs(dest / label, {
            "eunha2.exe": BINARY,
            "params.dat": PARENT / label / "params.dat",
            "khout.000001.dat": parent,
            "initial.dat": reference,
            **{name: DOCS / name for name in HELPERS},
        })
        if cases[label]["params.dat"] != "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43":
            raise ValueError(f"{label}: physics parameter file changed")
    manifest = {
        "run": str(dest), "start_time": 6.501, "target_time": 10,
        "binary_sha256": binary_sha, "mask_sha256": root_inputs["disk_active_mask.bin"],
        "parent_sha256": {label: accepted[label][0]["endpoint_sha256"] for label in accepted},
        "initial_sha256": {label: accepted[label][1]["initial_sha256"] for label in accepted},
        "root_inputs": root_inputs, "case_inputs": cases,
        "scope": "Two independent native weighted Kepler continuations from accepted t=6.501. "
                 "The weighted physics remains experimental pending local conservation and long-time gates.",
    }
    with (dest / "stage_manifest.json").open("x") as output:
        json.dump(manifest, output, indent=2)
        output.write("\n")
    print(json.dumps({"run": str(dest), "binary_sha256": binary_sha,
                      "parent_sha256": manifest["parent_sha256"]}, indent=2))


if __name__ == "__main__":
    main()

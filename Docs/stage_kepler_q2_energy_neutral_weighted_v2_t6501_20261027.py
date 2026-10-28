"""Stage v2 weighted native pair from both accepted t=6.5 checkpoints."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
BUILD = ROOT / "gpu_energy_neutral_native_build_v2"
BINARY = BUILD / "eunha2_energy_neutral_v2.exe"
PARAMS_SHA = "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"
RUNNER = "kepler_q2_energy_neutral_weighted_t6501_a10_20261027.slurm"
HELPERS = (
    "audit_kepler_q2_weighted_short_replay_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "laguerre_nonlinear_volume_probe_t6_20261026.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)
CASES = {
    "aligned": ("to65_near32_a10_v1", "KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json"),
    "jitter": ("jitter005_to65_a10_v1", "KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json"),
}


def copy_inputs(dest, files):
    dest.mkdir()
    checksums = {}
    for name, source in files.items():
        shutil.copy2(source, dest / name)
        checksums[name] = digest(dest / name)
    with (dest / "inputs.sha256").open("x") as stream:
        for name, value in sorted(checksums.items()):
            stream.write(f"{value}  {name}\n")
    return checksums


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("destination", type=Path)
    args = ap.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)
    binary_record = (BUILD / "binary.sha256").read_text().split()
    if len(binary_record) != 2 or Path(binary_record[1]) != BINARY:
        raise ValueError("Unexpected v2 build receipt")
    binary_sha = digest(BINARY)
    if binary_sha != binary_record[0]:
        raise ValueError("v2 native executable changed")
    patch = json.loads((BUILD / "patch_receipt.json").read_text())
    if patch["parent_source_sha256"] != "5ebc8d575bbb8bb361069923384b450682f950b36625c82c7396e50cb1c33f67":
        raise ValueError("Unexpected source patch parent")
    mask_receipt = json.loads((DOCS / "KEPLER_Q2_BOUNDED_NATIVE_WDOT_PROBE_INPUT_T65_2026-10-27.json").read_text())
    mask = Path(mask_receipt["mask"])
    if mask.stat().st_size != 55351 or digest(mask) != mask_receipt["mask_sha256"]:
        raise ValueError("Disk activity mask changed")
    parents = {}
    for name, (subdir, receipt_name) in CASES.items():
        receipt = json.loads((DOCS / receipt_name).read_text())
        source = ROOT / subdir
        endpoint = Path(receipt["endpoint"])
        if (receipt["time"] != 6.5 or len(receipt["checks"]) != 17
                or not all(receipt["checks"].values())
                or endpoint.parent != source
                or digest(endpoint) != receipt["endpoint_sha256"]
                or abs(snapshot_time(endpoint, 55351) - 6.5) > 1e-10
                or digest(source / "params.dat") != PARAMS_SHA):
            raise ValueError(f"{name}: parent gate or physics settings changed")
        parents[name] = (source, endpoint, receipt)
    dest.mkdir(parents=True)
    shutil.copy2(mask, dest / "disk_active_mask.bin")
    shutil.copy2(DOCS / RUNNER, dest / "job.slurm")
    with (dest / "inputs.sha256").open("x") as stream:
        for name in ("disk_active_mask.bin", "job.slurm"):
            stream.write(f"{digest(dest / name)}  {name}\n")
    staged = {}
    for name, (source, endpoint, _) in parents.items():
        staged[name] = copy_inputs(dest / name, {
            "eunha2.exe": BINARY,
            "params.dat": source / "params.dat",
            "khout.000001.dat": endpoint,
            **{helper: DOCS / helper for helper in HELPERS},
        })
    manifest = {
        "run": str(dest), "time_start": 6.5, "time_end": 6.501,
        "binary_sha256": binary_sha,
        "source_patch_receipt": patch,
        "mask_sha256": mask_receipt["mask_sha256"],
        "parent_sha256": {name: parents[name][2]["endpoint_sha256"] for name in CASES},
        "case_inputs": staged,
        "scope": "Experimental fixed-mass weighted RK4 one-step gate only; long Kepler remains unproven.",
    }
    with (dest / "stage_manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"run": str(dest), "binary_sha256": binary_sha,
                      "parent_sha256": manifest["parent_sha256"]}, indent=2))


if __name__ == "__main__":
    main()

"""Stage one-step W0 native FaceCSR Jacobi diagnostics from accepted states."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
BUILD = ROOT / "gpu_jacobi_native_build_v1/eunha2_wdot_probe.exe"
BINARY_SHA = "b731719d4455e7ec20fadb09e6c38d69e22391cf42f58d36e19ae108c80097ed"
PARAMS_SHA = "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"
RUNNER = "kepler_q2_native_gpu_wdot_probe_t6501_a10_20261027.slurm"
HELPERS = (
    "audit_kepler_bounded_short_replay_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)
CASES = {
    "aligned": ("to65_near32_a10_v1", "KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json"),
    "jitter": ("jitter005_to65_a10_v1", "KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("destination", type=Path)
    args = ap.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)
    if digest(BUILD) != BINARY_SHA:
        raise ValueError("Rebuilt native binary changed")
    mask_manifest = json.loads((DOCS / "KEPLER_Q2_BOUNDED_NATIVE_WDOT_PROBE_INPUT_T65_2026-10-27.json").read_text())
    mask = Path(mask_manifest["mask"])
    reference = Path(mask_manifest["rates"])
    if (mask.stat().st_size != 55351 or digest(mask) != mask_manifest["mask_sha256"]
            or digest(reference) != mask_manifest["rates_sha256"]):
        raise ValueError("Disk mask or CPU reference changed")
    sources = {}
    for label, (subdir, receipt_name) in CASES.items():
        receipt = json.loads((DOCS / receipt_name).read_text())
        if (receipt["time"] != 6.5 or len(receipt["checks"]) != 17
                or not all(receipt["checks"].values())
                or digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]
                or abs(snapshot_time(Path(receipt["endpoint"]), 55351)-6.5)>1e-10):
            raise ValueError(f"{label}: accepted parent gate changed")
        source = ROOT / subdir
        if (Path(receipt["endpoint"]).parent != source
                or digest(source / "params.dat") != PARAMS_SHA):
            raise ValueError(f"{label}: parent path/physics settings changed")
        sources[label] = (source, receipt)
    dest.mkdir(parents=True)
    root_inputs = {
        "disk_active_mask.bin": mask,
        "cpu_reference.npz": reference,
        "job.slurm": DOCS / RUNNER,
    }
    for target, source in root_inputs.items():
        shutil.copy2(source, dest / target)
    with (dest / "inputs.sha256").open("x") as stream:
        for target in sorted(root_inputs):
            stream.write(f"{digest(dest / target)}  {target}\n")
    staged = {}
    for label, (source, receipt) in sources.items():
        directory = dest / label
        directory.mkdir()
        inputs = {
            "eunha2.exe": BUILD,
            "params.dat": source / "params.dat",
            "khout.000001.dat": Path(receipt["endpoint"]),
            **{name: DOCS / name for name in HELPERS},
        }
        staged[label] = {}
        for target, original in inputs.items():
            shutil.copy2(original, directory / target)
            staged[label][target] = digest(directory / target)
        with (directory / "inputs.sha256").open("x") as stream:
            for target, sha in sorted(staged[label].items()):
                stream.write(f"{sha}  {target}\n")
    manifest = {
        "run": str(dest), "time_start": 6.5, "time_end": 6.501,
        "binary_sha256": BINARY_SHA,
        "mask_sha256": mask_manifest["mask_sha256"],
        "reference_sha256": mask_manifest["rates_sha256"],
        "accepted_parent_sha256": {label: sources[label][1]["endpoint_sha256"] for label in CASES},
        "inputs": staged,
        "scope": "Diagnostic-only GPU Jacobi16 from the native generated FaceCSR; W0 forces and RK4 state are not intentionally altered.",
    }
    with (dest / "stage_manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"run": str(dest), "binary_sha256": BINARY_SHA,
                      "mask_sha256": mask_manifest["mask_sha256"]}, indent=2))


if __name__ == "__main__":
    main()

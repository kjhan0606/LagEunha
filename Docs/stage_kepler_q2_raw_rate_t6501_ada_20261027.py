"""Replay one accepted q2 step using the isolated raw Jacobi weight rate."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
SOURCE = ROOT / "gpu_energy_neutral_capacity64_v5_fresh_localwork_t6501_ada_v1"
BUILD = ROOT / "gpu_energy_neutral_raw_rate_capacity64_build_v11b"
DEST = ROOT / "gpu_energy_neutral_raw_rate_v11b_t6501_ada_v1"
RECEIPT = DOCS / "KEPLER_FRESH_WEIGHT_LOCAL_T6501_2026-10-27/audit_t6501.json"


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    accepted = json.loads(RECEIPT.read_text())
    source_manifest = json.loads((SOURCE / "stage_manifest.json").read_text())
    build_manifest = json.loads((BUILD / "stage_manifest.json").read_text())
    if (accepted["time"] != 6.501 or accepted["parent_time"] != 6.5
            or len(accepted["checks"]) != 21 or not all(accepted["checks"].values())
            or accepted["binary_sha256"] != source_manifest["binary_sha256"]
            or accepted["parent_sha256"] != source_manifest["parent_endpoint_sha256"]
            or build_manifest["base"] != str(ROOT / "gpu_energy_neutral_capacity64_build_v5")):
        raise ValueError("Accepted v5 lineage does not match")
    for key, sha_key in (("initial", "initial_sha256"),
                         ("parent", "parent_sha256"),
                         ("endpoint", "endpoint_sha256")):
        if digest(Path(accepted[key])) != accepted[sha_key]:
            raise ValueError(f"Accepted {key} snapshot changed")
    binary = BUILD / "eunha2_raw_rate_v11.exe"
    binary_sha = digest(binary)
    if binary_sha != (BUILD / "binary.sha256").read_text().split()[0]:
        raise ValueError("Build binary hash mismatch")
    names = tuple(name for name in source_manifest["inputs"]
                  if name not in ("eunha2.exe", "job_ada.sh"))
    for name in names:
        if digest(SOURCE / name) != source_manifest["inputs"][name]:
            raise ValueError(f"Changed source input: {name}")
    runner = (SOURCE / "job_ada.sh").read_text()
    edits = {
        "export GFS_GPU_VOLW2=1 GFS_GPU_VOLW2_ITERS=64":
            "export GFS_GPU_VOLW2=1 GFS_GPU_VOLW2_ITERS=64 GFS_GPU_VOLW2_PROJECTION=0",
        "--expected-binary-sha256 " + accepted["binary_sha256"]:
            "--expected-binary-sha256 " + binary_sha,
    }
    for old, new in edits.items():
        if runner.count(old) != 1:
            raise ValueError(f"Unexpected runner anchor: {old}")
        runner = runner.replace(old, new)
    DEST.mkdir()
    hashes = {}
    for name in names:
        shutil.copy2(SOURCE / name, DEST / name)
        hashes[name] = digest(DEST / name)
        if hashes[name] != source_manifest["inputs"][name]:
            raise ValueError(f"Copy changed: {name}")
    shutil.copy2(binary, DEST / "eunha2.exe")
    hashes["eunha2.exe"] = digest(DEST / "eunha2.exe")
    if hashes["eunha2.exe"] != binary_sha:
        raise ValueError("Copied binary changed")
    (DEST / "job_ada.sh").write_text(runner)
    subprocess.run(["bash", "-n", str(DEST / "job_ada.sh")], check=True)
    hashes["job_ada.sh"] = digest(DEST / "job_ada.sh")
    with (DEST / "inputs.sha256").open("x") as out:
        for name, sha in sorted(hashes.items()):
            out.write(f"{sha}  {name}\n")
    with (DEST / "stage_manifest.json").open("x") as out:
        json.dump({"destination": str(DEST), "source_receipt": str(RECEIPT),
                   "source_receipt_sha256": digest(RECEIPT),
                   "parent_endpoint_sha256": accepted["parent_sha256"],
                   "baseline_endpoint_sha256": accepted["endpoint_sha256"],
                   "binary_build": str(BUILD), "binary_sha256": binary_sha,
                   "time_start": 6.5, "time_end": 6.501, "inputs": hashes,
                   "scope": "One native Ada step with raw Jacobi weight rate; "
                            "global pressure-work projection disabled only here."},
                  out, indent=2)
        out.write("\n")
    print(DEST)


if __name__ == "__main__":
    main()

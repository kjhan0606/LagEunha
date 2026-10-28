"""Test whether raw-rate volume control prevents q2 pressure noise from t=.825."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
REFERENCE = ROOT / "gpu_energy_neutral_capacity64_v5_recover_to65_a10_v1"
BUILD = ROOT / "gpu_energy_neutral_raw_rate_capacity64_build_v11b"
DEST = ROOT / "gpu_energy_neutral_raw_rate_v11b_from0825_to65_h200_v1"


def replace_once(source, old, new):
    if source.count(old) != 1:
        raise ValueError(f"Unexpected native runner anchor: {old!r}")
    return source.replace(old, new)


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    reference = json.loads((REFERENCE / "stage_manifest.json").read_text())
    target_audit = json.loads((REFERENCE / "audit_t65.json").read_text())
    parent_dir = Path(reference["parent"])
    parent_receipt = json.loads((parent_dir / reference["parent_receipt"]).read_text())
    if (reference["start"] != .825 or reference["stop"] != 6.5
            or reference["binary_sha256"] != parent_receipt["binary_sha256"]
            or parent_receipt["time"] != .825
            or target_audit["time"] != 6.5
            or target_audit["parent_sha256"] != parent_receipt["endpoint_sha256"]
            or not all(parent_receipt["checks"].values())
            or not all(target_audit["checks"].values())
            or digest(Path(parent_receipt["endpoint"])) != parent_receipt["endpoint_sha256"]):
        raise ValueError("Early weighted reference lineage is not accepted")
    for name, sha in reference["inputs"].items():
        if digest(REFERENCE / name) != sha:
            raise ValueError(f"Changed reference input: {name}")
    binary = BUILD / "eunha2_raw_rate_v11.exe"
    new_sha = digest(binary)
    if new_sha != (BUILD / "binary.sha256").read_text().split()[0]:
        raise ValueError("Raw-rate binary hash mismatch")
    old_sha = reference["binary_sha256"]
    runner = (REFERENCE / "job.slurm").read_text()
    for old, new in (
        ("LagQ2W64v5Chain2", "LagQ2RawFrom0825"),
        ("#SBATCH --partition=a10", "#SBATCH --partition=h200"),
        ("#SBATCH --gres=gpu:A10:1", "#SBATCH --gres=gpu:H200:1"),
        ("#SBATCH --time=03:00:00", "#SBATCH --time=06:00:00"),
        ("#SBATCH --nice=200", "#SBATCH --nice=10000"),
        ("export GFS_GPU_VOLW2=1 GFS_GPU_VOLW2_ITERS=64",
         "export GFS_GPU_VOLW2=1 GFS_GPU_VOLW2_ITERS=64 GFS_GPU_VOLW2_PROJECTION=0"),
        ("--expected-binary-sha256 " + old_sha,
         "--expected-binary-sha256 " + new_sha),
    ):
        runner = replace_once(runner, old, new)
    DEST.mkdir()
    hashes = {}
    for name, sha in reference["inputs"].items():
        if name in ("job.slurm", "eunha2.exe"):
            continue
        shutil.copy2(REFERENCE / name, DEST / name)
        hashes[name] = digest(DEST / name)
        if hashes[name] != sha:
            raise ValueError(f"Copied input changed: {name}")
    shutil.copy2(binary, DEST / "eunha2.exe")
    hashes["eunha2.exe"] = digest(DEST / "eunha2.exe")
    if hashes["eunha2.exe"] != new_sha:
        raise ValueError("Copied binary changed")
    (DEST / "job.slurm").write_text(runner)
    subprocess.run(["bash", "-n", str(DEST / "job.slurm")], check=True)
    hashes["job.slurm"] = digest(DEST / "job.slurm")
    with (DEST / "inputs.sha256").open("x") as out:
        for name, sha in sorted(hashes.items()):
            out.write(f"{sha}  {name}\n")
    with (DEST / "stage_manifest.json").open("x") as out:
        json.dump({"destination": str(DEST), "projected_reference": str(REFERENCE),
                   "projected_reference_receipt_sha256": digest(REFERENCE / "audit_t65.json"),
                   "accepted_parent_receipt": str(parent_dir / reference["parent_receipt"]),
                   "accepted_parent_endpoint_sha256": parent_receipt["endpoint_sha256"],
                   "binary_sha256": new_sha, "start": .825, "stop": 6.5,
                   "inputs": hashes,
                   "scope": "Independent H200 native raw-rate continuation from accepted "
                            "t=.825; checks whether first-ring pressure and torque "
                            "noise are prevented before the t=6.5 inherited state."},
                  out, indent=2)
        out.write("\n")
    print(DEST)


if __name__ == "__main__":
    main()

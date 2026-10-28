"""Stage an independent raw-rate H200 t=6.5->10 native continuation."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PROJECTED = ROOT / "gpu_energy_neutral_capacity64_v5_fresh_to10_a10_v1"
BUILD = ROOT / "gpu_energy_neutral_raw_rate_capacity64_build_v11b"
DEST = ROOT / "gpu_energy_neutral_raw_rate_v11b_t65_to10_h200_v1"
ONE_STEP = ROOT / "gpu_energy_neutral_raw_rate_v11b_t6501_ada_v1/audit_t6501.json"


def replace_once(source, old, new):
    if source.count(old) != 1:
        raise ValueError(f"Unexpected H200 runner anchor: {old!r}")
    return source.replace(old, new)


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    original_manifest = json.loads((PROJECTED / "stage_manifest.json").read_text())
    projected_receipt = json.loads((PROJECTED / "audit_t10.json").read_text())
    parent_dir = Path(original_manifest["parent"])
    parent_receipt = json.loads((parent_dir / original_manifest["parent_receipt"]).read_text())
    short = json.loads(ONE_STEP.read_text())
    if (original_manifest["start"] != 6.5 or original_manifest["stop"] != 10.
            or original_manifest["binary_sha256"] != parent_receipt["binary_sha256"]
            or parent_receipt["time"] != 6.5
            or parent_receipt["endpoint_sha256"] != short["parent_sha256"]
            or parent_receipt["initial_sha256"] != short["initial_sha256"]
            or not all(parent_receipt["checks"].values())
            or projected_receipt["time"] != 10.
            or not all(projected_receipt["checks"].values())):
        raise ValueError("Projected long-run lineage is not accepted")
    if digest(Path(parent_receipt["endpoint"])) != parent_receipt["endpoint_sha256"]:
        raise ValueError("Parent checkpoint changed")
    for name, sha in original_manifest["inputs"].items():
        if digest(PROJECTED / name) != sha:
            raise ValueError(f"Projected staged input changed: {name}")
    new_binary = BUILD / "eunha2_raw_rate_v11.exe"
    new_sha = digest(new_binary)
    if new_sha != (BUILD / "binary.sha256").read_text().split()[0]:
        raise ValueError("Raw-rate binary hash mismatch")
    old_sha = original_manifest["binary_sha256"]
    runner = (PROJECTED / "job.slurm").read_text()
    for old, new in (
        ("LagQ2W64v5Chain1", "LagQ2RawRate10"),
        ("#SBATCH --partition=a10", "#SBATCH --partition=h200"),
        ("#SBATCH --gres=gpu:A10:1", "#SBATCH --gres=gpu:H200:1"),
        ("#SBATCH --nice=200", "#SBATCH --nice=10000"),
        ("export GFS_GPU_VOLW2=1 GFS_GPU_VOLW2_ITERS=64",
         "export GFS_GPU_VOLW2=1 GFS_GPU_VOLW2_ITERS=64 GFS_GPU_VOLW2_PROJECTION=0"),
        ("--expected-binary-sha256 " + old_sha,
         "--expected-binary-sha256 " + new_sha),
    ):
        runner = replace_once(runner, old, new)
    DEST.mkdir()
    hashes = {}
    for name, sha in original_manifest["inputs"].items():
        if name in ("job.slurm", "eunha2.exe"):
            continue
        shutil.copy2(PROJECTED / name, DEST / name)
        hashes[name] = digest(DEST / name)
        if hashes[name] != sha:
            raise ValueError(f"Copy changed: {name}")
    shutil.copy2(new_binary, DEST / "eunha2.exe")
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
        json.dump({"destination": str(DEST), "projected_run": str(PROJECTED),
                   "projected_receipt_sha256": digest(PROJECTED / "audit_t10.json"),
                   "parent_receipt": str(parent_dir / original_manifest["parent_receipt"]),
                   "parent_endpoint_sha256": parent_receipt["endpoint_sha256"],
                   "short_raw_receipt_sha256": digest(ONE_STEP),
                   "binary_sha256": new_sha, "start": 6.5, "stop": 10.,
                   "inputs": hashes,
                   "scope": "One-rank H200 native raw-rate continuation from accepted "
                            "t=6.5, compared with completed projected t=10; "
                            "no t=0 replay and no modification of existing jobs."},
                  out, indent=2)
        out.write("\n")
    print(DEST)


if __name__ == "__main__":
    main()

"""Stage the v11b default-projection control from the same t=6.5 state."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
SOURCE = ROOT / "gpu_energy_neutral_raw_rate_v11b_t6501_ada_v1"
DEST = ROOT / "gpu_energy_neutral_projected_v11b_t6501_ada_v1"


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    source = json.loads((SOURCE / "stage_manifest.json").read_text())
    receipt = json.loads((SOURCE / "audit_t6501.json").read_text())
    if (receipt["parent_sha256"] != source["parent_endpoint_sha256"]
            or receipt["binary_sha256"] != source["binary_sha256"]
            or len(receipt["checks"]) != 21 or not all(receipt["checks"].values())):
        raise ValueError("Raw pilot lineage or gate failed")
    for name, sha in source["inputs"].items():
        if digest(SOURCE / name) != sha:
            raise ValueError(f"Changed input: {name}")
    runner = (SOURCE / "job_ada.sh").read_text()
    old = "GFS_GPU_VOLW2_PROJECTION=0"
    if runner.count(old) != 1:
        raise ValueError("Expected one raw-rate option")
    runner = runner.replace(old, "GFS_GPU_VOLW2_PROJECTION=1")
    DEST.mkdir()
    hashes = {}
    for name, sha in source["inputs"].items():
        if name == "job_ada.sh":
            continue
        shutil.copy2(SOURCE / name, DEST / name)
        hashes[name] = digest(DEST / name)
        if hashes[name] != sha:
            raise ValueError(f"Copy changed: {name}")
    (DEST / "job_ada.sh").write_text(runner)
    subprocess.run(["bash", "-n", str(DEST / "job_ada.sh")], check=True)
    hashes["job_ada.sh"] = digest(DEST / "job_ada.sh")
    with (DEST / "inputs.sha256").open("x") as out:
        for name, sha in sorted(hashes.items()):
            out.write(f"{sha}  {name}\n")
    with (DEST / "stage_manifest.json").open("x") as out:
        json.dump({"destination": str(DEST), "parent_raw_pilot": str(SOURCE),
                   "parent_endpoint_sha256": source["parent_endpoint_sha256"],
                   "binary_sha256": source["binary_sha256"],
                   "inputs": hashes,
                   "scope": "Same v11b binary and t=6.5 parent, restored default "
                            "pressure-work projection as a bitwise control."},
                  out, indent=2)
        out.write("\n")
    print(DEST)


if __name__ == "__main__":
    main()

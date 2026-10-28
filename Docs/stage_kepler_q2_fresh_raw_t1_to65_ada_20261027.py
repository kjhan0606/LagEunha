"""Continue accepted same-t0 fresh raw q2 checkpoint from t=1 to 6.5."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
SOURCE = ROOT / "gpu_energy_neutral_raw_rate_v11b_t1_to65_ada_v1"
FRESH = ROOT / "gpu_energy_neutral_v11b_fresh_raw_t001_to1_a10_v1"
DEST = ROOT / "gpu_energy_neutral_raw_rate_v11b_fresh_t1_to65_ada_v1"


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    source_manifest = json.loads((SOURCE / "stage_manifest.json").read_text())
    source_audit = json.loads((SOURCE / "audit_t65.json").read_text())
    parent = json.loads((FRESH / "audit_t1.json").read_text())
    if (source_audit["time"] != 6.5 or len(source_audit["checks"]) != 21
            or not all(source_audit["checks"].values())
            or parent["time"] != 1. or parent["parent_time"] != .001
            or len(parent["checks"]) != 21 or not all(parent["checks"].values())
            or parent["binary_sha256"] != source_manifest["binary_sha256"]
            or parent["initial_sha256"] != source_audit["initial_sha256"]):
        raise ValueError("Accepted fresh parent or prior Ada template changed")
    for key, sha_key in (("initial", "initial_sha256"),
                         ("parent", "parent_sha256"),
                         ("endpoint", "endpoint_sha256")):
        if digest(Path(parent[key])) != parent[sha_key]:
            raise ValueError(f"Changed fresh {key} snapshot")
    for name, sha in source_manifest["inputs"].items():
        if digest(SOURCE / name) != sha:
            raise ValueError(f"Changed source input: {name}")
    runner = (SOURCE / "job_ada.sh").read_text()
    if (runner.count("GFS_GPU_VOLW2_PROJECTION=0") != 1
            or runner.count("HYDRO_TSTOP=6.5 EUNHA_DUMP_DT=0.825") != 1
            or runner.count("--parent-time 1.0 --time 6.5") != 1
            or runner.count("--output audit_t65.json") != 1):
        raise ValueError("Unexpected raw Ada runner")
    DEST.mkdir()
    hashes = {}
    for name, sha in source_manifest["inputs"].items():
        if name in ("job_ada.sh", "khout.000001.dat"):
            continue
        shutil.copy2(SOURCE / name, DEST / name)
        hashes[name] = digest(DEST / name)
        if hashes[name] != sha:
            raise ValueError(f"Changed copied input: {name}")
    if hashes["initial.dat"] != parent["initial_sha256"]:
        raise ValueError("Fresh raw t=0 reference differs from source")
    shutil.copy2(Path(parent["endpoint"]), DEST / "khout.000001.dat")
    hashes["khout.000001.dat"] = digest(DEST / "khout.000001.dat")
    if hashes["khout.000001.dat"] != parent["endpoint_sha256"]:
        raise ValueError("Copied fresh t=1 checkpoint changed")
    (DEST / "job_ada.sh").write_text(runner)
    subprocess.run(["bash", "-n", str(DEST / "job_ada.sh")], check=True)
    hashes["job_ada.sh"] = digest(DEST / "job_ada.sh")
    with (DEST / "inputs.sha256").open("x") as output:
        for name, sha in sorted(hashes.items()):
            output.write(f"{sha}  {name}\n")
    with (DEST / "stage_manifest.json").open("x") as output:
        json.dump({"destination": str(DEST),
                   "accepted_parent_receipt": str(FRESH / "audit_t1.json"),
                   "accepted_parent_endpoint_sha256": parent["endpoint_sha256"],
                   "binary_sha256": parent["binary_sha256"],
                   "start": 1., "stop": 6.5, "inputs": hashes,
                   "scope": "Same-binary fresh raw trajectory from accepted t=0 → .001 "
                            "→ 1 checkpoint; no projected .825 history."},
                  output, indent=2)
        output.write("\n")
    print(DEST)


if __name__ == "__main__":
    main()

"""Continue the accepted raw-rate Ada q2 trajectory from t=6.5 to 10."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
SOURCE = ROOT / "gpu_energy_neutral_raw_rate_v11b_t1_to65_ada_v1"
DEST = ROOT / "gpu_energy_neutral_raw_rate_v11b_t65_to10_ada_v1"


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    manifest = json.loads((SOURCE / "stage_manifest.json").read_text())
    accepted = json.loads((SOURCE / "audit_t65.json").read_text())
    if (accepted["time"] != 6.5 or accepted["parent_time"] != 1.
            or len(accepted["checks"]) != 21 or not all(accepted["checks"].values())
            or accepted["binary_sha256"] != manifest["binary_sha256"]
            or accepted["parent_sha256"] != manifest["accepted_parent_endpoint_sha256"]):
        raise ValueError("Accepted raw t=6.5 parent changed")
    for key, sha_key in (("initial", "initial_sha256"),
                         ("parent", "parent_sha256"),
                         ("endpoint", "endpoint_sha256")):
        if digest(Path(accepted[key])) != accepted[sha_key]:
            raise ValueError(f"Changed {key} snapshot")
    for name, sha in manifest["inputs"].items():
        if digest(SOURCE / name) != sha:
            raise ValueError(f"Changed source input: {name}")
    runner = (SOURCE / "job_ada.sh").read_text()
    changes = {
        "HYDRO_TSTOP=6.5 EUNHA_DUMP_DT=0.825":
            "HYDRO_TSTOP=10.0 EUNHA_DUMP_DT=0.875",
        "--parent-time 1.0 --time 6.5":
            "--parent-time 6.5 --time 10.0",
        "--output audit_t65.json": "--output audit_t10.json",
    }
    for old, new in changes.items():
        if runner.count(old) != 1:
            raise ValueError(f"Unexpected runner anchor: {old}")
        runner = runner.replace(old, new)
    if "GFS_GPU_VOLW2_PROJECTION=0" not in runner:
        raise ValueError("Raw-rate option missing")
    DEST.mkdir()
    hashes = {}
    for name, sha in manifest["inputs"].items():
        if name in ("job_ada.sh", "khout.000001.dat"):
            continue
        shutil.copy2(SOURCE / name, DEST / name)
        hashes[name] = digest(DEST / name)
        if hashes[name] != sha:
            raise ValueError(f"Changed copy: {name}")
    shutil.copy2(Path(accepted["endpoint"]), DEST / "khout.000001.dat")
    hashes["khout.000001.dat"] = digest(DEST / "khout.000001.dat")
    if hashes["khout.000001.dat"] != accepted["endpoint_sha256"]:
        raise ValueError("Copied accepted t=6.5 checkpoint changed")
    (DEST / "job_ada.sh").write_text(runner)
    subprocess.run(["bash", "-n", str(DEST / "job_ada.sh")], check=True)
    hashes["job_ada.sh"] = digest(DEST / "job_ada.sh")
    with (DEST / "inputs.sha256").open("x") as output:
        for name, sha in sorted(hashes.items()):
            output.write(f"{sha}  {name}\n")
    with (DEST / "stage_manifest.json").open("x") as output:
        json.dump({"destination": str(DEST),
                   "accepted_parent_receipt": str(SOURCE / "audit_t65.json"),
                   "accepted_parent_endpoint_sha256": accepted["endpoint_sha256"],
                   "binary_sha256": accepted["binary_sha256"],
                   "start": 6.5, "stop": 10., "inputs": hashes,
                   "scope": "One Ada native raw-rate continuation from accepted "
                            "t=6.5 endpoint; no fresh t=0 replay."},
                  output, indent=2)
        output.write("\n")
    print(DEST)


if __name__ == "__main__":
    main()

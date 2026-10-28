"""Continue accepted raw-rate t=1 Ada checkpoint to t=6.5 without a t0 replay."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
SOURCE = ROOT / "gpu_energy_neutral_raw_rate_v11b_t0825_to1_ada_v1"
DEST = ROOT / "gpu_energy_neutral_raw_rate_v11b_t1_to65_ada_v1"


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    source_manifest = json.loads((SOURCE / "stage_manifest.json").read_text())
    accepted = json.loads((SOURCE / "audit_t1.json").read_text())
    if (source_manifest["mode"] != "raw" or accepted["time"] != 1.
            or accepted["parent_time"] != .825
            or accepted["binary_sha256"] != source_manifest["binary_sha256"]
            or accepted["parent_sha256"] != source_manifest["parent_endpoint_sha256"]
            or len(accepted["checks"]) != 21 or not all(accepted["checks"].values())):
        raise ValueError("Accepted t=1 raw parent is not intact")
    for key, sha_key in (("initial", "initial_sha256"),
                         ("parent", "parent_sha256"),
                         ("endpoint", "endpoint_sha256")):
        if digest(Path(accepted[key])) != accepted[sha_key]:
            raise ValueError(f"Changed {key} snapshot")
    for name, sha in source_manifest["inputs"].items():
        if digest(SOURCE / name) != sha:
            raise ValueError(f"Changed source input: {name}")
    runner = (SOURCE / "job_ada.sh").read_text()
    changes = {
        "HYDRO_TSTOP=1.0 EUNHA_DUMP_DT=0.175":
            "HYDRO_TSTOP=6.5 EUNHA_DUMP_DT=0.825",
        "--parent-time 0.825 --time 1.0":
            "--parent-time 1.0 --time 6.5",
        "--output audit_t1.json": "--output audit_t65.json",
    }
    for old, new in changes.items():
        if runner.count(old) != 1:
            raise ValueError(f"Unexpected runner anchor: {old}")
        runner = runner.replace(old, new)
    if "GFS_GPU_VOLW2_PROJECTION=0" not in runner:
        raise ValueError("Raw-rate option missing")
    DEST.mkdir()
    hashes = {}
    for name, sha in source_manifest["inputs"].items():
        if name in ("job_ada.sh", "khout.000001.dat"):
            continue
        shutil.copy2(SOURCE / name, DEST / name)
        hashes[name] = digest(DEST / name)
        if hashes[name] != sha:
            raise ValueError(f"Copied input changed: {name}")
    shutil.copy2(Path(accepted["endpoint"]), DEST / "khout.000001.dat")
    hashes["khout.000001.dat"] = digest(DEST / "khout.000001.dat")
    if hashes["khout.000001.dat"] != accepted["endpoint_sha256"]:
        raise ValueError("Copied t=1 checkpoint changed")
    (DEST / "job_ada.sh").write_text(runner)
    subprocess.run(["bash", "-n", str(DEST / "job_ada.sh")], check=True)
    hashes["job_ada.sh"] = digest(DEST / "job_ada.sh")
    with (DEST / "inputs.sha256").open("x") as out:
        for name, sha in sorted(hashes.items()):
            out.write(f"{sha}  {name}\n")
    with (DEST / "stage_manifest.json").open("x") as out:
        json.dump({"destination": str(DEST), "accepted_parent_receipt": str(SOURCE / "audit_t1.json"),
                   "accepted_parent_endpoint_sha256": accepted["endpoint_sha256"],
                   "binary_sha256": accepted["binary_sha256"],
                   "start": 1., "stop": 6.5, "inputs": hashes,
                   "scope": "One Ada GPU raw-weight native continuation from "
                            "accepted t=1 endpoint; no t=0 restart. "
                            "Diagnostic alternative while H200 QOS is saturated."},
                  out, indent=2)
        out.write("\n")
    print(DEST)


if __name__ == "__main__":
    main()

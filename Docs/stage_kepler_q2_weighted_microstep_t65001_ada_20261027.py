"""Stage a guarded 0.0001 native weighted microstep from accepted t=6.5."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
SOURCE = ROOT/"gpu_energy_neutral_capacity64_v5_fresh_localwork_t6501_ada_v1"
DEST = ROOT/"gpu_energy_neutral_capacity64_v5_microstep_t65001_ada_v1"
RECEIPT = DOCS/"KEPLER_FRESH_WEIGHT_LOCAL_T6501_2026-10-27/audit_t6501.json"


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    accepted = json.loads(RECEIPT.read_text())
    source_manifest = json.loads((SOURCE/"stage_manifest.json").read_text())
    if (len(accepted["checks"]) != 21 or not all(accepted["checks"].values())
            or accepted["parent_time"] != 6.5
            or accepted["parent_sha256"] != source_manifest["parent_endpoint_sha256"]
            or accepted["binary_sha256"] != source_manifest["binary_sha256"]):
        raise ValueError("Accepted weighted parent is not intact")
    for key, sha_key in (("initial", "initial_sha256"),
                         ("parent", "parent_sha256"),
                         ("endpoint", "endpoint_sha256")):
        if digest(Path(accepted[key])) != accepted[sha_key]:
            raise ValueError(f"Changed {key} snapshot")
    for name, sha in source_manifest["inputs"].items():
        if digest(SOURCE/name) != sha:
            raise ValueError(f"Changed source input: {name}")
    runner = (SOURCE/"job_ada.sh").read_text()
    changes = {
        "HYDRO_TSTOP=6.501 EUNHA_DUMP_DT=0.001":
            "HYDRO_TSTOP=6.5001 EUNHA_DUMP_DT=0.0001",
        "--parent-time 6.5 --time 6.501":
            "--parent-time 6.5 --time 6.5001",
        "--output audit_t6501.json": "--output audit_t65001.json",
    }
    for old, new in changes.items():
        if runner.count(old) != 1:
            raise ValueError(f"Unexpected runner anchor: {old}")
        runner = runner.replace(old, new)
    DEST.mkdir()
    hashes = {}
    for name, sha in source_manifest["inputs"].items():
        if name == "job_ada.sh":
            continue
        shutil.copy2(SOURCE/name, DEST/name)
        hashes[name] = digest(DEST/name)
        if hashes[name] != sha:
            raise ValueError(f"Copy changed: {name}")
    (DEST/"job_ada.sh").write_text(runner)
    subprocess.run(["bash", "-n", str(DEST/"job_ada.sh")], check=True)
    hashes["job_ada.sh"] = digest(DEST/"job_ada.sh")
    with (DEST/"inputs.sha256").open("x") as stream:
        for name, sha in sorted(hashes.items()):
            stream.write(f"{sha}  {name}\n")
    report = {"destination": str(DEST), "source_receipt": str(RECEIPT),
              "source_receipt_sha256": digest(RECEIPT),
              "parent_endpoint_sha256": accepted["parent_sha256"],
              "source_runner_sha256": source_manifest["inputs"]["job_ada.sh"],
              "time_start": 6.5, "time_end": 6.5001,
              "inputs": hashes,
              "scope": "Same accepted v5 weighted checkpoint and binary; one shorter native RK4 step on Ada."}
    with (DEST/"stage_manifest.json").open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(DEST)


if __name__ == "__main__":
    main()

"""Stage two same-binary Ada continuations from the accepted t=6.5 state."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
CASES = {
    "projected": (ROOT / "gpu_energy_neutral_projected_v11b_t6501_ada_v1",
                  ROOT / "gpu_energy_neutral_projected_v11b_t65_to66_ada_v1"),
    "raw": (ROOT / "gpu_energy_neutral_raw_rate_v11b_t6501_ada_v1",
            ROOT / "gpu_energy_neutral_raw_rate_v11b_t65_to66_ada_v1"),
}


def main():
    receipts = {name: json.loads((source / "audit_t6501.json").read_text())
                for name, (source, _) in CASES.items()}
    if (receipts["projected"]["parent_sha256"] != receipts["raw"]["parent_sha256"]
            or receipts["projected"]["initial_sha256"] != receipts["raw"]["initial_sha256"]
            or receipts["projected"]["binary_sha256"] != receipts["raw"]["binary_sha256"]):
        raise ValueError("Same-binary pair lineage differs")
    for name, (source, dest) in CASES.items():
        if dest.exists():
            raise FileExistsError(dest)
        manifest = json.loads((source / "stage_manifest.json").read_text())
        receipt = receipts[name]
        if (receipt["time"] != 6.501 or receipt["parent_time"] != 6.5
                or receipt["parent_sha256"] != manifest["parent_endpoint_sha256"]
                or receipt["binary_sha256"] != manifest["binary_sha256"]
                or len(receipt["checks"]) != 21 or not all(receipt["checks"].values())):
            raise ValueError(f"Unaccepted short {name} replay")
        for filename, sha in manifest["inputs"].items():
            if digest(source / filename) != sha:
                raise ValueError(f"Changed {name} input: {filename}")
        runner = (source / "job_ada.sh").read_text()
        for old, new in (
            ("HYDRO_TSTOP=6.501 EUNHA_DUMP_DT=0.001",
             "HYDRO_TSTOP=6.6 EUNHA_DUMP_DT=0.1"),
            ("--parent-time 6.5 --time 6.501", "--parent-time 6.5 --time 6.6"),
            ("--output audit_t6501.json", "--output audit_t66.json"),
        ):
            if runner.count(old) != 1:
                raise ValueError(f"Runner anchor changed: {old}")
            runner = runner.replace(old, new)
        dest.mkdir()
        hashes = {}
        for filename, sha in manifest["inputs"].items():
            if filename == "job_ada.sh":
                continue
            shutil.copy2(source / filename, dest / filename)
            hashes[filename] = digest(dest / filename)
            if hashes[filename] != sha:
                raise ValueError(f"Copied {name} input changed: {filename}")
        (dest / "job_ada.sh").write_text(runner)
        subprocess.run(["bash", "-n", str(dest / "job_ada.sh")], check=True)
        hashes["job_ada.sh"] = digest(dest / "job_ada.sh")
        with (dest / "inputs.sha256").open("x") as out:
            for filename, sha in sorted(hashes.items()):
                out.write(f"{sha}  {filename}\n")
        with (dest / "stage_manifest.json").open("x") as out:
            json.dump({"destination": str(dest), "mode": name,
                       "source_receipt": str(source / "audit_t6501.json"),
                       "source_receipt_sha256": digest(source / "audit_t6501.json"),
                       "parent_endpoint_sha256": receipt["parent_sha256"],
                       "binary_sha256": receipt["binary_sha256"],
                       "time_start": 6.5, "time_end": 6.6,
                       "inputs": hashes,
                       "scope": "Same v11b binary and accepted parent; isolated "
                                "100-step Ada comparison of weight projection."},
                      out, indent=2)
            out.write("\n")
        print(dest)


if __name__ == "__main__":
    main()

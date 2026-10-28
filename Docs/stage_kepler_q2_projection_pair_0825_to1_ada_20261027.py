"""Stage same-binary Ada projection pair from the accepted t=.825 checkpoint."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PARENT = ROOT / "gpu_energy_neutral_capacity64_v5_recover_to0825_a10_v1/audit_t0825.json"
CASES = {
    "projected": (ROOT / "gpu_energy_neutral_projected_v11b_t6501_ada_v1",
                  ROOT / "gpu_energy_neutral_projected_v11b_t0825_to1_ada_v1"),
    "raw": (ROOT / "gpu_energy_neutral_raw_rate_v11b_t6501_ada_v1",
            ROOT / "gpu_energy_neutral_raw_rate_v11b_t0825_to1_ada_v1"),
}


def main():
    accepted = json.loads(PARENT.read_text())
    if (accepted["time"] != .825 or len(accepted["checks"]) != 21
            or not all(accepted["checks"].values())
            or digest(Path(accepted["endpoint"])) != accepted["endpoint_sha256"]):
        raise ValueError("t=.825 parent is not accepted")
    for name, (source, dest) in CASES.items():
        if dest.exists():
            raise FileExistsError(dest)
        manifest = json.loads((source / "stage_manifest.json").read_text())
        receipt = json.loads((source / "audit_t6501.json").read_text())
        if (receipt["time"] != 6.501 or not all(receipt["checks"].values())
                or receipt["initial_sha256"] != accepted["initial_sha256"]
                or receipt["binary_sha256"] != manifest["binary_sha256"]
                or manifest["inputs"]["disk_active_mask.bin"] !=
                   digest(source / "disk_active_mask.bin")):
            raise ValueError(f"Invalid {name} v11b source lineage")
        for filename, sha in manifest["inputs"].items():
            if digest(source / filename) != sha:
                raise ValueError(f"Changed {name} input: {filename}")
        runner = (source / "job_ada.sh").read_text()
        for old, new in (
            ("HYDRO_TSTOP=6.501 EUNHA_DUMP_DT=0.001",
             "HYDRO_TSTOP=1.0 EUNHA_DUMP_DT=0.175"),
            ("--parent-time 6.5 --time 6.501", "--parent-time 0.825 --time 1.0"),
            ("--output audit_t6501.json", "--output audit_t1.json"),
        ):
            if runner.count(old) != 1:
                raise ValueError(f"Runner anchor changed: {old}")
            runner = runner.replace(old, new)
        dest.mkdir()
        hashes = {}
        for filename, sha in manifest["inputs"].items():
            if filename in ("job_ada.sh", "khout.000001.dat"):
                continue
            shutil.copy2(source / filename, dest / filename)
            hashes[filename] = digest(dest / filename)
            if hashes[filename] != sha:
                raise ValueError(f"Copied {name} input changed: {filename}")
        shutil.copy2(Path(accepted["endpoint"]), dest / "khout.000001.dat")
        hashes["khout.000001.dat"] = digest(dest / "khout.000001.dat")
        if hashes["khout.000001.dat"] != accepted["endpoint_sha256"]:
            raise ValueError("Copied t=.825 parent changed")
        (dest / "job_ada.sh").write_text(runner)
        subprocess.run(["bash", "-n", str(dest / "job_ada.sh")], check=True)
        hashes["job_ada.sh"] = digest(dest / "job_ada.sh")
        with (dest / "inputs.sha256").open("x") as out:
            for filename, sha in sorted(hashes.items()):
                out.write(f"{sha}  {filename}\n")
        with (dest / "stage_manifest.json").open("x") as out:
            json.dump({"destination": str(dest), "mode": name,
                       "accepted_parent_receipt": str(PARENT),
                       "parent_endpoint_sha256": accepted["endpoint_sha256"],
                       "binary_sha256": manifest["binary_sha256"],
                       "start": .825, "stop": 1., "inputs": hashes,
                       "scope": "Same Ada GPU and v11b binary; independent "
                                "early raw/projected moving-weight control."},
                      out, indent=2)
            out.write("\n")
        print(dest)


if __name__ == "__main__":
    main()

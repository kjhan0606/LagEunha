"""Continue both accepted fresh t=.001 q2 branches to t=1 on native A10."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PAIR = DOCS / "KEPLER_Q2_FRESH_PROJECTION_PAIR_T001_2026-10-27.json"
INITIAL_SHA = "8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b"
BINARY_SHA = "5d10fd4c17ac67bfafcc461fda03937bc645159758d1fff1ab73bd27ecf54f72"


def main():
    pair = json.loads(PAIR.read_text())
    if (pair["time"] != .001 or len(pair["initial_state_fields_exact"]) != 16
            or any(pair["modes"][m]["audit_passed"] != 18 for m in ("projected", "raw"))):
        raise ValueError("Fresh raw/projected pair is not accepted")
    template = (DOCS / "kepler_q2_fresh_projection_to1_20261027.slurm").read_text()
    if template.count("__BINARY_SHA__") != 2:
        raise ValueError("Unexpected continuation job template")
    for mode in ("projected", "raw"):
        if (ROOT / f"gpu_energy_neutral_v11b_fresh_{mode}_t001_to1_a10_v1").exists():
            raise FileExistsError(mode)

    for mode, projection in (("projected", 1), ("raw", 0)):
        parent = ROOT / f"gpu_energy_neutral_v11b_fresh_{mode}_t001_a10_v1"
        dest = ROOT / f"gpu_energy_neutral_v11b_fresh_{mode}_t001_to1_a10_v1"
        receipt = json.loads((parent / "audit_t001.json").read_text())
        if (receipt["time"] != .001 or len(receipt["checks"]) != 18
                or not all(receipt["checks"].values())
                or receipt["binary_sha256"] != BINARY_SHA
                or receipt["reference_sha256"] != INITIAL_SHA
                or receipt["endpoint_sha256"] != pair["modes"][mode]["endpoint_sha256"]
                or digest(parent / "khout.000001.dat") != receipt["endpoint_sha256"]):
            raise ValueError(f"Unaccepted {mode} fresh parent")
        sources = {
            "eunha2.exe": parent / "eunha2.exe",
            "params.dat": parent / "params.dat",
            "disk_active_mask.bin": parent / "disk_active_mask.bin",
            "initial.dat": parent / "accepted_initial.dat",
        }
        for name in (
            "audit_kepler_q2_weighted_equilibrium_20261027.py",
            "audit_kepler_bounded_native_evolution_20261027.py",
            "audit_kepler_sparse_native_ballistic_20261027.py",
            "audit_kepler_sparse_area_t0825_20261027.py",
            "audit_kepler_sparse_state_consistency_20261027.py",
            "laguerre_nonlinear_volume_probe_t6_20261026.py",
            "prepare_kepler_sparse_area_restart_20261027.py",
        ):
            sources[name] = DOCS / name
        if digest(sources["eunha2.exe"]) != BINARY_SHA or digest(sources["initial.dat"]) != INITIAL_SHA:
            raise ValueError("Binary or t=0 reference changed")
        dest.mkdir()
        hashes = {}
        for name, source in sources.items():
            shutil.copy2(source, dest / name)
            hashes[name] = digest(dest / name)
            if hashes[name] != digest(source):
                raise ValueError(f"Changed copy: {name}")
        runner = (template.replace("__MODE__", mode[:3])
                  .replace("__PARENT_DIR__", str(parent))
                  .replace("__INITIAL_SHA__", INITIAL_SHA)
                  .replace("__BINARY_SHA__", BINARY_SHA)
                  .replace("__PROJECTION__", str(projection)))
        if "__" in runner:
            raise ValueError("Unresolved continuation job placeholder")
        (dest / "job.slurm").write_text(runner)
        subprocess.run(["bash", "-n", str(dest / "job.slurm")], check=True)
        hashes["job.slurm"] = digest(dest / "job.slurm")
        with (dest / "inputs.sha256").open("x") as output:
            for name, sha in sorted(hashes.items()):
                output.write(f"{sha}  {name}\n")
        with (dest / "stage_manifest.json").open("x") as output:
            json.dump({"run": str(dest), "mode": mode, "projection": projection,
                       "accepted_parent_receipt": str(parent / "audit_t001.json"),
                       "accepted_parent_endpoint_sha256": receipt["endpoint_sha256"],
                       "initial_sha256": INITIAL_SHA, "binary_sha256": BINARY_SHA,
                       "start": .001, "stop": 1., "inputs": hashes,
                       "scope": "Same-binary fresh t=0 lineage continuation from "
                                "each accepted t=.001 raw/projected checkpoint."},
                      output, indent=2)
            output.write("\n")
        print(dest)


if __name__ == "__main__":
    main()

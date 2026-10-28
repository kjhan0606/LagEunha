"""Stage an exact-time, same-IC t=10 comparison and density rendering."""
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
RUN = ROOT / "fresh_compare_density_t10_cpu_v1"
W0_RECEIPT = ROOT / "jitter005_to10_a10_v1" / "audit_t10.json"
WEIGHTED_RECEIPT = (
    ROOT / "gpu_energy_neutral_capacity64_v5_fresh_to10_a10_v1" / "audit_t10.json"
)
FILES = (
    "kepler_fresh_weighted_compare_density_t10_20261027.slurm",
    "audit_kepler_bounded_entropy_sign_20261027.py",
    "audit_kepler_weighted_entropy_sign_20261027.py",
    "audit_kepler_weighted_vs_w0_fresh_20261027.py",
    "plot_kepler_q2_weighted_vs_w0_density_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)
EXPECTED_W0_FAILURES = {"disk_ballistic_residual"}


def main():
    if RUN.exists():
        raise FileExistsError(RUN)
    w0 = json.loads(W0_RECEIPT.read_text())
    weighted = json.loads(WEIGHTED_RECEIPT.read_text())
    w0_failed = {key for key, passed in w0["checks"].items() if not passed}
    weighted_failed = {key for key, passed in weighted["checks"].items() if not passed}
    if (w0["initial_sha256"] != weighted["initial_sha256"]
            or w0["time"] != 10.0 or weighted["time"] != 10.0
            or w0_failed != EXPECTED_W0_FAILURES or weighted_failed
            or len(weighted["checks"]) != 21):
        raise ValueError("Unexpected same-IC t=10 receipts or gate results")

    subprocess.run(["bash", "-n", str(DOCS / FILES[0])], check=True)
    RUN.mkdir()
    hashes = {}
    for name in FILES:
        source = DOCS / name
        destination = RUN / name
        shutil.copy2(source, destination)
        hashes[name] = digest(destination)
        if hashes[name] != digest(source):
            raise ValueError(f"Changed staged file: {name}")
    for source, name in ((W0_RECEIPT, "w0_receipt.json"),
                         (WEIGHTED_RECEIPT, "weighted_receipt.json")):
        destination = RUN / name
        shutil.copy2(source, destination)
        hashes[name] = digest(destination)
        if hashes[name] != digest(source):
            raise ValueError(f"Changed receipt: {source}")
    with (RUN / "inputs.sha256").open("x") as output:
        for name, sha256 in sorted(hashes.items()):
            output.write(f"{sha256}  {name}\n")
    with (RUN / "stage_manifest.json").open("x") as output:
        json.dump({"scope": "Same initial snapshot, exact t=10 endpoints; W0's one failed ballistic-residual gate is retained and labeled.",
                   "initial_sha256": w0["initial_sha256"],
                   "w0_endpoint_sha256": w0["endpoint_sha256"],
                   "weighted_endpoint_sha256": weighted["endpoint_sha256"],
                   "input_sha256": hashes}, output, indent=2)
        output.write("\n")
    print(RUN)


if __name__ == "__main__":
    main()

"""Stage the true from-early-time weighted versus W0 t10 density image."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DEST = ROOT / "fresh_weighted_vs_w0_density_t10_cpu_v2"
W0 = ROOT / "jitter005_to10_a10_v1/audit_t10.json"
WEIGHTED = ROOT / "gpu_energy_neutral_capacity64_v5_fresh_to10_a10_v1/audit_t10.json"
NAMES = ("kepler_q2_fresh_weighted_vs_w0_density_t10_cpu_20261027.slurm",
         "plot_kepler_q2_weighted_vs_w0_density_20261027.py",
         "audit_kepler_sparse_area_t0825_20261027.py",
         "prepare_kepler_sparse_area_restart_20261027.py")


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    a, b = json.loads(W0.read_text()), json.loads(WEIGHTED.read_text())
    if (a["time"] != 10. or b["time"] != 10.
            or a["initial_sha256"] != b["initial_sha256"]
            or {k for k, value in a["checks"].items() if not value}
               != {"disk_ballistic_residual"}
            or not all(b["checks"].values())
            or digest(Path(a["endpoint"])) != a["endpoint_sha256"]
            or digest(Path(b["endpoint"])) != b["endpoint_sha256"]):
        raise ValueError("Fresh weighted/W0 t10 lineage or gates differ")
    runner = DOCS / NAMES[0]
    subprocess.run(["bash", "-n", str(runner)], check=True)
    DEST.mkdir()
    hashes = {}
    for name in NAMES:
        shutil.copy2(DOCS / name, DEST / name)
        hashes[name] = digest(DEST / name)
    (DEST / "job.slurm").symlink_to(NAMES[0])
    with (DEST / "inputs.sha256").open("x") as out:
        for name, sha in sorted(hashes.items()):
            out.write(f"{sha}  {name}\n")
    with (DEST / "stage_manifest.json").open("x") as out:
        json.dump({"destination": str(DEST), "w0_receipt": str(W0),
                   "weighted_receipt": str(WEIGHTED),
                   "w0_endpoint_sha256": a["endpoint_sha256"],
                   "weighted_endpoint_sha256": b["endpoint_sha256"],
                   "initial_sha256": a["initial_sha256"],
                   "inputs": hashes,
                   "scope": "Actual surface-density and relative-density fields "
                            "from same-IC native t10 W0 and fresh weighted runs."},
                  out, indent=2)
        out.write("\n")
    print(DEST)


if __name__ == "__main__":
    main()

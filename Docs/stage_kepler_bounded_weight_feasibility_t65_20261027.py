"""Stage two independent CPU-only static Laguerre volume feasibility audits."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
DEST = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/weight_feasibility_t65_cpu_v1")
NAMES = (
    "audit_kepler_bounded_weight_feasibility_t65_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "laguerre_newton_volume_t6_20261026.py",
    "laguerre_nonlinear_volume_probe_t6_20261026.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
    "kepler_q2_bounded_weight_feasibility_t65_cpu_20261027.slurm",
)
RECEIPTS = {"aligned_receipt.json": "KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json",
            "jitter_receipt.json": "KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json"}


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    subprocess.run(["bash", "-n", str(DOCS / NAMES[-1])], check=True)
    sources = {name: DOCS / name for name in NAMES}
    sources.update({target: DOCS / source for target, source in RECEIPTS.items()})
    for target in RECEIPTS:
        receipt = json.loads(sources[target].read_text())
        if (receipt["time"] != 6.5 or len(receipt["checks"]) != 17
                or not all(receipt["checks"].values())
                or digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]):
            raise ValueError(f"Invalid accepted receipt: {target}")
    DEST.mkdir()
    for target, source in sources.items():
        shutil.copy2(source, DEST / target)
    with (DEST / "source.sha256").open("x") as stream:
        for target in sorted(sources):
            stream.write(f"{digest(DEST / target)}  {target}\n")
    manifest = {"run": str(DEST), "time": 6.5,
                "work": "Two independent actual-position static Laguerre disk-volume solves; no hydro evolution",
                "files": {target: digest(DEST / target) for target in sorted(sources)}}
    (DEST / "staging.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

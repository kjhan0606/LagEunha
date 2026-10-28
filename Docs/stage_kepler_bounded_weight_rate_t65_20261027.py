"""Stage independent one-stage static Laguerre weight-rate probes."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
SOURCE = ROOT / "weight_feasibility_t65_cpu_v1"
DEST = ROOT / "weight_rate_t65_cpu_v1"
SCRIPTS = (
    "audit_kepler_bounded_weight_rate_t65_20261027.py",
    "audit_kepler_bounded_weight_feasibility_t65_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "laguerre_newton_volume_t6_20261026.py",
    "laguerre_nonlinear_volume_probe_t6_20261026.py",
    "laguerre_weight_rate_actual_t625_20261026.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
    "kepler_q2_bounded_weight_rate_t65_cpu_20261027.slurm",
)


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    subprocess.run(["bash", "-n", str(DOCS / SCRIPTS[-1])], check=True)
    files = {name: DOCS / name for name in SCRIPTS}
    for variant in ("aligned", "jitter"):
        result = json.loads((SOURCE / f"{variant}_weight.json").read_text())
        if (not result["converged_to_1e_minus_7"]
                or result["native_face_guard"]["faces_over_0p8"]):
            raise ValueError(f"Static weighted geometry not accepted: {variant}")
        files[f"{variant}_receipt.json"] = SOURCE / f"{variant}_receipt.json"
        files[f"{variant}_weight.npz"] = SOURCE / f"{variant}_weight.npz"
    DEST.mkdir()
    for target, source in files.items():
        shutil.copy2(source, DEST / target)
    with (DEST / "source.sha256").open("x") as stream:
        for target in sorted(files):
            stream.write(f"{digest(DEST / target)}  {target}\n")
    manifest = {"run": str(DEST), "time": 6.5,
                "source_static_solve": str(SOURCE),
                "work": "One geometry-only weight-rate solve and 0.00193983-time-unit Euler prediction per IC",
                "files": {target: digest(DEST / target) for target in sorted(files)}}
    (DEST / "staging.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

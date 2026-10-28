"""Stage the disk-face material-speed least-squares weight-rate comparison."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
SOURCE = ROOT / "weighted_face_slip_t65_cpu_v1"
DEST = ROOT / "material_rate_t65_cpu_v1"
SCRIPTS = (
    "audit_kepler_bounded_material_rate_t65_20261027.py",
    "audit_kepler_bounded_weight_feasibility_t65_20261027.py",
    "audit_kepler_bounded_weighted_face_slip_t65_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "laguerre_nonlinear_volume_probe_t6_20261026.py",
    "laguerre_weight_rate_actual_t625_20261026.py",
    "laguerre_newton_volume_t6_20261026.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
    "kepler_q2_bounded_material_rate_t65_cpu_20261027.slurm",
)


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    subprocess.run(["bash", "-n", str(DOCS / SCRIPTS[-1])], check=True)
    files = {name: DOCS / name for name in SCRIPTS}
    for variant in ("aligned", "jitter"):
        receipt = SOURCE / f"{variant}_receipt.json"
        weight = SOURCE / f"{variant}_weight.npz"
        rate = SOURCE / f"{variant}_rate.npz"
        prior = json.loads((SOURCE / f"{variant}_slip.json").read_text())
        if (prior["time"] != 6.5 or digest(weight) != prior["static_weights_sha256"]
                or digest(rate) != prior["weight_rates_sha256"]
                or json.loads(receipt.read_text())["endpoint_sha256"] != prior["endpoint_sha256"]):
            raise ValueError(f"Incompatible material-rate inputs: {variant}")
        files[f"{variant}_receipt.json"] = receipt
        files[f"{variant}_weight.npz"] = weight
        files[f"{variant}_rate.npz"] = rate
    DEST.mkdir()
    for target, source in files.items():
        shutil.copy2(source, DEST / target)
    with (DEST / "source.sha256").open("x") as stream:
        for target in sorted(files):
            stream.write(f"{digest(DEST / target)}  {target}\n")
    manifest = {"run": str(DEST), "time": 6.5,
                "work": "Constrained-volume rate versus disk-face material-speed least-squares rate",
                "files": {target: digest(DEST / target) for target in sorted(files)}}
    (DEST / "staging.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

"""Stage independent material-face speed comparisons at bounded t=6.5."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
STATIC = ROOT / "weight_feasibility_t65_cpu_v1"
RATE = ROOT / "weight_rate_t65_cpu_v1"
DEST = ROOT / "weighted_face_slip_t65_cpu_v1"
SCRIPTS = (
    "audit_kepler_bounded_weighted_face_slip_t65_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "laguerre_nonlinear_volume_probe_t6_20261026.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
    "kepler_q2_bounded_weighted_face_slip_t65_cpu_20261027.slurm",
)


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    subprocess.run(["bash", "-n", str(DOCS / SCRIPTS[-1])], check=True)
    files = {name: DOCS / name for name in SCRIPTS}
    for variant in ("aligned", "jitter"):
        result = json.loads((RATE / f"{variant}_rate.json").read_text())
        weight = STATIC / f"{variant}_weight.npz"
        receipt = STATIC / f"{variant}_receipt.json"
        if (result["time"] != 6.5 or digest(weight) != result["static_weights_sha256"]
                or json.loads(receipt.read_text())["endpoint_sha256"] != result["endpoint_sha256"]):
            raise ValueError(f"Incompatible rate and static inputs: {variant}")
        files[f"{variant}_receipt.json"] = receipt
        files[f"{variant}_weight.npz"] = weight
        files[f"{variant}_rate.npz"] = RATE / f"{variant}_rate.npz"
    DEST.mkdir()
    for target, source in files.items():
        shutil.copy2(source, DEST / target)
    with (DEST / "source.sha256").open("x") as stream:
        for target in sorted(files):
            stream.write(f"{digest(DEST / target)}  {target}\n")
    manifest = {"run": str(DEST), "time": 6.5,
                "work": "Analytic Kepler material velocity versus Voronoi/Laguerre face normal velocity",
                "files": {target: digest(DEST / target) for target in sorted(files)}}
    (DEST / "staging.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

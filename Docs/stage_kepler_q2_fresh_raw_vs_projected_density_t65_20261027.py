"""Stage accepted fresh-raw versus projected exact-t6.5 density image."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS=Path(__file__).resolve().parent
ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DEST=ROOT/"fresh_raw_vs_projected_density_t65_inner_cpu_v1"
BASE=ROOT/"gpu_energy_neutral_capacity64_v5_recover_to65_a10_v1/audit_t65.json"
RAW=ROOT/"gpu_energy_neutral_raw_rate_v11b_fresh_t1_to65_ada_v1/audit_t65.json"
NAMES=("kepler_q2_fresh_raw_vs_projected_density_t65_cpu_20261027.slurm",
       "plot_kepler_q2_projected_vs_raw_density_20261027.py",
       "plot_kepler_q2_weighted_vs_w0_density_20261027.py",
       "audit_kepler_sparse_area_t0825_20261027.py",
       "prepare_kepler_sparse_area_restart_20261027.py")


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    a,b=json.loads(BASE.read_text()),json.loads(RAW.read_text())
    if (a["time"]!=6.5 or b["time"]!=6.5
            or a["initial_sha256"]!=b["initial_sha256"]
            or len(a["checks"])!=21 or len(b["checks"])!=21
            or not all(a["checks"].values()) or not all(b["checks"].values())
            or digest(Path(a["endpoint"]))!=a["endpoint_sha256"]
            or digest(Path(b["endpoint"]))!=b["endpoint_sha256"]):
        raise ValueError("Exact native t=6.5 pair not accepted")
    runner=DOCS/NAMES[0]
    subprocess.run(["bash","-n",str(runner)],check=True)
    DEST.mkdir()
    hashes={}
    for name in NAMES:
        shutil.copy2(DOCS/name,DEST/name)
        hashes[name]=digest(DEST/name)
    (DEST/"job.slurm").symlink_to(NAMES[0])
    with (DEST/"inputs.sha256").open("x") as output:
        for name,sha in sorted(hashes.items()):
            output.write(f"{sha}  {name}\n")
    with (DEST/"stage_manifest.json").open("x") as output:
        json.dump({"destination":str(DEST),"projected_receipt":str(BASE),
                   "raw_receipt":str(RAW),
                   "projected_endpoint_sha256":a["endpoint_sha256"],
                   "raw_endpoint_sha256":b["endpoint_sha256"],
                   "initial_sha256":a["initial_sha256"],"inputs":hashes,
                   "scope":"Same-t0, exact-t6.5 native power-cell surface-density "
                           "and residual fields with common color scales and inner zoom."},
                  output,indent=2)
        output.write("\n")
    print(DEST)


if __name__=="__main__":
    main()

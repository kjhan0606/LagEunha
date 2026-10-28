"""Prestage a dependent coherent-ring W0 continuation to exact t=6.5."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS=Path(__file__).resolve().parent
PARENT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/ringcoherent005_t001_a10_v1")
RUNNER="kepler_q2_ring_coherent_to65_a10_20261027.slurm"
SCRIPTS={
    "job.slurm":RUNNER,
    "prepare_kepler_q2_ring_from_preflight_20261027.py":"prepare_kepler_q2_ring_from_preflight_20261027.py",
    "audit_kepler_bounded_native_evolution_20261027.py":"audit_kepler_bounded_native_evolution_20261027.py",
    "audit_kepler_sparse_native_ballistic_20261027.py":"audit_kepler_sparse_native_ballistic_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py":"audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py":"audit_kepler_sparse_state_consistency_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py":"prepare_kepler_sparse_area_restart_20261027.py",
}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("destination",type=Path)
    ap.add_argument("--parent-job",type=int,required=True)
    args=ap.parse_args()
    dest=args.destination.resolve()
    if dest.exists(): raise FileExistsError(dest)
    if args.parent_job<=0: raise ValueError("Invalid parent job")
    subprocess.run(["bash","-n",str(DOCS/RUNNER)],check=True)
    parent_manifest=json.loads((PARENT/"stage_manifest.json").read_text())
    if parent_manifest["run"]!=str(PARENT): raise ValueError("Wrong parent")
    dest.mkdir(parents=True)
    hashes={}
    for target,source in SCRIPTS.items():
        shutil.copy2(DOCS/source,dest/target)
        hashes[target]=digest(dest/target)
    with (dest/"scripts.sha256").open("x") as stream:
        for name,value in sorted(hashes.items()):
            stream.write(f"{value}  {name}\n")
    manifest={"run":str(dest),"parent":str(PARENT),"parent_job":args.parent_job,
              "dependency":f"afterok:{args.parent_job}","scripts":hashes,
              "scope":"Run only after coherent-ring preflight passes; native W0 t=6.5 IC-sensitivity control."}
    with (dest/"prestage_manifest.json").open("x") as stream:
        json.dump(manifest,stream,indent=2)
        stream.write("\n")
    print(json.dumps(manifest,indent=2))


if __name__=="__main__": main()

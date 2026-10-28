"""Stage an IC-only coherent-ring control with the accepted W0 native binary."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS=Path(__file__).resolve().parent
ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
SOURCE=ROOT/"jitter005_t001_a10_v1"
PREFIX="KEPLER_Q2_BOUNDED_RINGCOHERENT005_R275_N1024_2026-10-27"
RUNNER="kepler_q2_ring_coherent_preflight_a10_20261027.slurm"
HELPERS=("audit_kepler_bounded_native_preflight_20261027.py",
         "audit_kepler_sparse_area_t0825_20261027.py",
         "audit_kepler_sparse_state_consistency_20261027.py",
         "prepare_kepler_sparse_area_restart_20261027.py")
BINARY_SHA="3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e"
PARAMS_SHA="2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("destination",type=Path)
    args=ap.parse_args()
    dest=args.destination.resolve()
    if dest.exists(): raise FileExistsError(dest)
    subprocess.run(["bash","-n",str(DOCS/RUNNER)],check=True)
    accepted=json.loads((SOURCE/"audit_t001.json").read_text())
    if accepted["n"]!=55351 or not all(accepted["checks"].values()):
        raise ValueError("Reference W0 preflight is not accepted")
    if digest(SOURCE/"eunha2.exe")!=BINARY_SHA or digest(SOURCE/"params.dat")!=PARAMS_SHA:
        raise ValueError("Accepted native binary/params changed")
    ic=json.loads((DOCS/(PREFIX+".json")).read_text())
    glass=DOCS/(PREFIX+".glass.bin")
    areas=DOCS/(PREFIX+".areas.bin")
    if (ic["sites"]!=55351 or ic["disk_sites"]!=49038 or
            ic["phase_mode"]!="one shared angle offset per ring" or
            digest(glass)!=ic["glass_sha256"] or digest(areas)!=ic["areas_sha256"]):
        raise ValueError("Coherent-ring IC changed")
    sources={"eunha2.exe":SOURCE/"eunha2.exe",
             "params.dat":SOURCE/"params.dat",
             "glass.bin":glass,"areas.bin":areas,
             "ic_manifest.json":DOCS/(PREFIX+".json"),
             "job.slurm":DOCS/RUNNER,
             **{name:DOCS/name for name in HELPERS}}
    dest.mkdir(parents=True)
    checksums={}
    for name,source in sources.items():
        shutil.copy2(source,dest/name)
        checksums[name]=digest(dest/name)
    with (dest/"inputs.sha256").open("x") as stream:
        for name,value in sorted(checksums.items()):
            stream.write(f"{value}  {name}\n")
    with (dest/"stage_manifest.json").open("x") as stream:
        json.dump({"run":str(dest),"reference":str(SOURCE),
                   "reference_audit_sha256":digest(SOURCE/"audit_t001.json"),
                   "ic_manifest_sha256":checksums["ic_manifest.json"],
                   "inputs":checksums,
                   "scope":"IC sensitivity preflight; coherent ring offsets versus aligned and independent site jitter."},
                  stream,indent=2)
        stream.write("\n")
    print(json.dumps({"run":str(dest),"glass_sha256":checksums["glass.bin"],
                      "areas_sha256":checksums["areas.bin"]},indent=2))


if __name__=="__main__": main()

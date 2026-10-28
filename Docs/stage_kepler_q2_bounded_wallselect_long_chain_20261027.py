"""Prestage conditional t=30→60→100 bounded W0 continuations."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DOCS=Path(__file__).resolve().parent
RUNNER="kepler_q2_bounded_wallselect_chain_20261027.slurm"
PREPARER="prepare_kepler_bounded_wallselect_chain_20261027.py"
SOURCE=ROOT/"wallselect_pair_t58_a10_v1/candidate/eunha2.exe"
SOURCE_SHA="87d60c1cb059e1710dd2749c84af5c98ba17fe077e6fe080b9c8b42d5fd40968"
DESTINATIONS=(ROOT/"wallselect_to60_a10_v1",ROOT/"wallselect_to100_a10_v1")


def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def main():
    if any(path.exists() for path in DESTINATIONS):
        raise FileExistsError("A target run directory already exists")
    if digest(SOURCE)!=SOURCE_SHA:
        raise ValueError("Validated selective-wall binary changed")
    subprocess.run(["bash","-n",str(DOCS/RUNNER)],check=True)
    for dest,start,target in zip(DESTINATIONS,(30,60),(60,100)):
        dest.mkdir()
        for name in (RUNNER,PREPARER):
            shutil.copy2(DOCS/name,dest/name)
            if digest(dest/name)!=digest(DOCS/name):
                raise ValueError(f"Copy changed: {name}")
        manifest={"destination":str(dest),"parent_time":start,"target_time":target,
                  "candidate_binary":str(SOURCE),"candidate_binary_sha256":SOURCE_SHA,
                  "prestaged_sha256":{name:digest(dest/name) for name in (RUNNER,PREPARER)},
                  "scope":"Accepted-parent-only bounded W0 continuation; no claim "
                          "of GIZMO-matched initial conditions or cell spin."}
        with (dest/"prestage_manifest.json").open("x") as output:
            json.dump(manifest,output,indent=2)
            output.write("\n")
        print(json.dumps({"destination":str(dest),"start":start,"target":target}))


if __name__=="__main__":
    main()

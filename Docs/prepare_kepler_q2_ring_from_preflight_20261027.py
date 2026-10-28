"""Copy an accepted coherent-ring preflight into a fresh W0 continuation."""

import json
from pathlib import Path
import shutil

from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PARENT=ROOT/"ringcoherent005_t001_a10_v1"
BINARY_SHA="3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e"
PARAMS_SHA="2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"
COPY={"eunha2.exe":"eunha2.exe","params.dat":"params.dat",
      "glass.bin":"glass.bin","areas.bin":"areas.bin",
      "initial.dat":"khout.000000.dat",
      "khout.000001.dat":"khout.000001.dat",
      "ic_manifest.json":"ic_manifest.json"}
SCRIPTS=(
    "job.slurm","prepare_kepler_q2_ring_from_preflight_20261027.py",
    "audit_kepler_bounded_native_evolution_20261027.py",
    "audit_kepler_sparse_native_ballistic_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)


def main():
    dest=Path.cwd().resolve()
    if any((dest/name).exists() for name in COPY) or (dest/"inputs.sha256").exists():
        raise FileExistsError("Continuation already contains physical inputs")
    audit=json.loads((PARENT/"audit_t001.json").read_text())
    if (audit["n"]!=55351 or abs(audit["t0"])>1e-12 or
            abs(audit["t1"]-.001)>1e-12 or not all(audit["checks"].values())):
        raise ValueError("Coherent-ring native preflight has not passed")
    for name in SCRIPTS:
        if not (dest/name).is_file(): raise FileNotFoundError(dest/name)
    if (digest(PARENT/"eunha2.exe")!=BINARY_SHA or
            digest(PARENT/"params.dat")!=PARAMS_SHA or
            digest(PARENT/"areas.bin")!=audit["area_sha256"] or
            digest(PARENT/"khout.000000.dat")!=audit["initial_sha256"] or
            digest(PARENT/"khout.000001.dat")!=audit["endpoint_sha256"] or
            abs(snapshot_time(PARENT/"khout.000001.dat",55351)-.001)>1e-12):
        raise ValueError("Preflight binary, sidecar, or state changed")
    ic=json.loads((PARENT/"ic_manifest.json").read_text())
    if (ic["phase_mode"]!="one shared angle offset per ring" or
            digest(PARENT/"glass.bin")!=ic["glass_sha256"] or
            digest(PARENT/"areas.bin")!=ic["areas_sha256"]):
        raise ValueError("Coherent-ring IC changed")
    for target,source in COPY.items():
        shutil.copy2(PARENT/source,dest/target)
    with (dest/"inputs.sha256").open("x") as stream:
        for name in sorted((*COPY,*SCRIPTS)):
            stream.write(f"{digest(dest/name)}  {name}\n")
    manifest={"parent":str(PARENT),"parent_audit_sha256":digest(PARENT/"audit_t001.json"),
              "accepted_parent_time":.001,"initial":str(dest/"initial.dat"),
              "initial_sha256":digest(dest/"initial.dat"),
              "parent_endpoint_sha256":audit["endpoint_sha256"],
              "binary_sha256":BINARY_SHA,"runner":"job.slurm"}
    with (dest/"restart_origin.json").open("x") as stream:
        json.dump(manifest,stream,indent=2)
        stream.write("\n")
    print(json.dumps(manifest,indent=2))


if __name__=="__main__": main()

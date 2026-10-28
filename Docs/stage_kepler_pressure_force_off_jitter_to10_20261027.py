"""Stage a t6.5→10 diagnostic from the same accepted jitter checkpoint."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


DOCS=Path(__file__).resolve().parent
ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PARENT=ROOT/"jitter005_to65_a10_v1"
SHORT=ROOT/"pressure_force_off_jitter_t65_short_a10_v1"
ORIGINAL=ROOT/"jitter005_to10_a10_v1"
DEST=ROOT/"pressure_force_off_jitter_to10_a10_v1"
BASE_SHA="3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e"
CANDIDATE_SHA="0a65c373613a5b31f8449b456387ebcf80d72136a81e46eec10c9692a8bed5e7"
RUNNER="kepler_pressure_force_off_jitter_to10_a10_20261027.slurm"
HELPERS=("audit_kepler_pressure_force_off_long_20261027.py",
         "audit_kepler_bounded_native_evolution_20261027.py",
         "audit_kepler_sparse_area_t0825_20261027.py",
         "audit_kepler_sparse_native_ballistic_20261027.py",
         "audit_kepler_sparse_state_consistency_20261027.py",
         "prepare_kepler_sparse_area_restart_20261027.py")


def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    parent=json.loads((PARENT/"audit_t65.json").read_text())
    original=json.loads((ORIGINAL/"audit_t10.json").read_text())
    short=json.loads((SHORT/"pressure_off_pair.json").read_text())
    if (parent["time"]!=6.5 or not all(parent["checks"].values())
            or original["time"]!=10
            or [k for k,v in original["checks"].items() if not v]!=["disk_ballistic_residual"]
            or not all(short["checks"].values())
            or abs(short["difference_acceleration_over_frozen_pressure_force_rms"]-1)>0.01
            or short["accepted_parent_sha256"]!=parent["endpoint_sha256"]
            or digest(Path(parent["endpoint"]))!=parent["endpoint_sha256"]
            or digest(Path(parent["initial"]))!=parent["initial_sha256"]
            or digest(SHORT/"eunha2.exe")!=CANDIDATE_SHA
            or digest(PARENT/"eunha2.exe")!=BASE_SHA):
        raise ValueError("Native parent or one-step causal control changed")
    subprocess.run(["bash","-n",str(DOCS/RUNNER)],check=True)
    DEST.mkdir()
    files={"eunha2.exe":SHORT/"eunha2.exe",
           "params.dat":PARENT/"params.dat",
           "glass.bin":PARENT/"glass.bin",
           "areas.bin":PARENT/"areas.bin",
           "khout.000001.dat":Path(parent["endpoint"]),
           RUNNER:DOCS/RUNNER,
           **{name:DOCS/name for name in HELPERS}}
    for name,source in files.items():
        shutil.copy2(source,DEST/name)
    manifest={"destination":str(DEST),"accepted_parent_receipt":str(PARENT/"audit_t65.json"),
              "accepted_parent_sha256":parent["endpoint_sha256"],
              "initial":str(Path(parent["initial"]).resolve()),
              "initial_sha256":parent["initial_sha256"],
              "one_step_pair":str(SHORT/"pressure_off_pair.json"),
              "one_step_pair_sha256":digest(SHORT/"pressure_off_pair.json"),
              "original_t10_receipt":str(ORIGINAL/"audit_t10.json"),
              "original_t10_receipt_sha256":digest(ORIGINAL/"audit_t10.json"),
              "diagnostic_binary_sha256":CANDIDATE_SHA,
              "source_files_sha256":{name:digest(DEST/name) for name in files},
              "scope":"Pressure-force-off causal control; geometric PdV retained; not physical hydro"}
    (DEST/"staging.json").write_text(json.dumps(manifest,indent=2)+"\n")
    with (DEST/"inputs.sha256").open("x") as stream:
        for name in sorted(files):
            stream.write(f"{digest(DEST/name)}  {name}\n")
        stream.write(f"{digest(DEST/'staging.json')}  staging.json\n")
    print(json.dumps({"destination":str(DEST),
                      "accepted_parent_sha256":parent["endpoint_sha256"],
                      "diagnostic_binary_sha256":CANDIDATE_SHA}))


if __name__=="__main__":
    main()

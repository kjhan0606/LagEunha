"""Stage one short native causal control from the accepted jitter t6.5 state."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


DOCS=Path(__file__).resolve().parent
ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PARENT=ROOT/"jitter005_to65_a10_v1"
DEST=ROOT/"pressure_force_off_jitter_t65_short_a10_v1"
BASE=ROOT/"force_inputs_t65_a10_v1/jitter"
CANDIDATE=Path("/tmp/Kepler_pressure_force_off_20261027/eunha2.exe")
CANDIDATE_SHA="0a65c373613a5b31f8449b456387ebcf80d72136a81e46eec10c9692a8bed5e7"
BASE_SHA="3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e"
RUNNER="kepler_pressure_force_off_jitter_short_a10_20261027.slurm"
HELPERS=("audit_kepler_bounded_short_replay_20261027.py",
         "audit_kepler_pressure_force_off_pair_20261027.py",
         "audit_kepler_sparse_area_t0825_20261027.py",
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
    r=json.loads((PARENT/"audit_t65.json").read_text())
    if (r["time"]!=6.5 or not all(r["checks"].values())
            or digest(Path(r["endpoint"]))!=r["endpoint_sha256"]
            or digest(Path(r["initial"]))!=r["initial_sha256"]
            or digest(PARENT/"eunha2.exe")!=BASE_SHA
            or digest(CANDIDATE)!=CANDIDATE_SHA):
        raise ValueError("Changed native parent or causal-control binary")
    b=json.loads((BASE/"audit_short.json").read_text())
    if not all(b["checks"].values()) or b["restart_sha256"]!=r["endpoint_sha256"]:
        raise ValueError("Original one-step reference is not accepted")
    subprocess.run(["bash","-n",str(DOCS/RUNNER)],check=True)
    DEST.mkdir()
    files={"eunha2.exe":CANDIDATE,
           "params.dat":PARENT/"params.dat",
           "glass.bin":PARENT/"glass.bin",
           "areas.bin":PARENT/"areas.bin",
           "khout.000001.dat":Path(r["endpoint"]),
           "pressure_force_t65.json":DOCS/"KEPLER_Q2_BOUNDED_PRESSURE_FORCE_T10_2026-10-27/jitter_t65.json",
           RUNNER:DOCS/RUNNER,
           **{name:DOCS/name for name in HELPERS}}
    for name,source in files.items():
        shutil.copy2(source,DEST/name)
    with (DEST/"inputs.sha256").open("x") as stream:
        for name in sorted(files):
            stream.write(f"{digest(DEST/name)}  {name}\n")
    manifest={"destination":str(DEST),"accepted_parent_receipt":str(PARENT/"audit_t65.json"),
              "accepted_parent_endpoint_sha256":r["endpoint_sha256"],
              "original_short_reference":str(BASE/"audit_short.json"),
              "original_short_endpoint_sha256":b["endpoint_sha256"],
              "original_binary_sha256":BASE_SHA,
              "diagnostic_binary_sha256":CANDIDATE_SHA,
              "files_sha256":{name:digest(DEST/name) for name in files},
              "scope":"One-step pressure-force-off intervention only; PdV and EOS retained; not physical hydro"}
    (DEST/"staging.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps({"destination":str(DEST),"parent_endpoint_sha256":r["endpoint_sha256"],
                      "diagnostic_binary_sha256":CANDIDATE_SHA}))


if __name__=="__main__":
    main()

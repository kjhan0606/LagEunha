"""Stage one read-only wall-face probe from accepted bounded t=.825 state."""

import hashlib
import json
from pathlib import Path
import shutil


DOCS = Path("/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs")
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PREVIOUS = ROOT/"to0825_near32_v1"
DEST = ROOT/"wall0825_probe_near32_a10_v1"
BINARY_SHA = "3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    receipt_path = PREVIOUS/"audit_t0825.json"
    audit = json.loads(receipt_path.read_text())
    if abs(audit["time"]-.825)>1e-10 or not all(audit["checks"].values()):
        raise ValueError("Accepted t=.825 gate missing")
    endpoint = Path(audit["endpoint"])
    if endpoint.parent.resolve()!=PREVIOUS.resolve() or digest(endpoint)!=audit["endpoint_sha256"]:
        raise ValueError("Accepted t=.825 checkpoint mismatch")
    if digest(PREVIOUS/"eunha2.exe")!=BINARY_SHA:
        raise ValueError("Accepted near32 binary mismatch")
    files = {
        "eunha2.exe": PREVIOUS/"eunha2.exe",
        "params.dat": PREVIOUS/"params.dat",
        "khout.000001.dat": endpoint,
        "job.slurm": DOCS/"kepler_q2_bounded_wall0825_probe_a10_20261027.slurm",
        "audit_kepler_fixedpoly_wall_faces_20261027.py":
            DOCS/"audit_kepler_fixedpoly_wall_faces_20261027.py",
        "audit_kepler_bounded_short_replay_20261027.py":
            DOCS/"audit_kepler_bounded_short_replay_20261027.py",
        "audit_kepler_sparse_area_t0825_20261027.py":
            DOCS/"audit_kepler_sparse_area_t0825_20261027.py",
        "audit_kepler_sparse_state_consistency_20261027.py":
            DOCS/"audit_kepler_sparse_state_consistency_20261027.py",
        "prepare_kepler_sparse_area_restart_20261027.py":
            DOCS/"prepare_kepler_sparse_area_restart_20261027.py",
    }
    hashes={name:digest(source) for name,source in files.items()}
    DEST.mkdir(exist_ok=False)
    for name,source in files.items():
        shutil.copy2(source,DEST/name)
        if digest(DEST/name)!=hashes[name]:
            raise IOError(f"Staged {name} changed")
    (DEST/"inputs.sha256").write_text("".join(
        f"{value}  {name}\n" for name,value in sorted(hashes.items())))
    stage={"run":str(DEST),"predecessor":str(PREVIOUS),
           "predecessor_audit_sha256":digest(receipt_path),
           "accepted_checkpoint_sha256":audit["endpoint_sha256"],
           "target_time":.826,"binary_sha256":BINARY_SHA,"input_sha256":hashes}
    (DEST/"staging.json").write_text(json.dumps(stage,indent=2)+"\n")
    print(json.dumps(stage,indent=2))


if __name__=="__main__":
    main()

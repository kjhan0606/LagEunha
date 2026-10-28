"""Stage separate CPU static weight solves for aligned and failed jitter t10."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS=Path(__file__).resolve().parent
DEST=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/weight_feasibility_t10_cpu_v1")
SOURCE=("audit_kepler_bounded_weight_feasibility_t65_20261027.py",
        "audit_kepler_sparse_area_t0825_20261027.py",
        "laguerre_newton_volume_t6_20261026.py",
        "laguerre_nonlinear_volume_probe_t6_20261026.py",
        "prepare_kepler_sparse_area_restart_20261027.py",
        "kepler_q2_bounded_weight_feasibility_t10_cpu_20261027.slurm")
RECEIPTS={"aligned_receipt.json":"KEPLER_Q2_BOUNDED_T10_DIAG_2026-10-27/aligned_audit_t10.json",
          "jitter_receipt.json":"KEPLER_Q2_BOUNDED_T10_DIAG_2026-10-27/jitter_failed_audit_t10.json"}
PARSER=Path("/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/audit_snapshot.py")


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    subprocess.run(["bash","-n",str(DOCS/SOURCE[-1])],check=True)
    for name,target in RECEIPTS.items():
        r=json.loads((DOCS/target).read_text())
        expected=set() if name.startswith("aligned") else {"disk_ballistic_residual"}
        failed={key for key,passed in r["checks"].items() if not passed}
        if (r["time"]!=10 or failed!=expected
                or digest(Path(r["initial"]))!=r["initial_sha256"]
                or digest(Path(r["endpoint"]))!=r["endpoint_sha256"]):
            raise ValueError(f"Native t10 receipt changed: {name}")
    DEST.mkdir()
    sources={name:DOCS/name for name in SOURCE}
    sources.update({name:DOCS/path for name,path in RECEIPTS.items()})
    for name,path in sources.items():
        shutil.copy2(path,DEST/name)
    with (DEST/"source.sha256").open("x") as stream:
        for name in sorted(sources):
            stream.write(f"{digest(DEST/name)}  {name}\n")
        stream.write(f"{digest(PARSER)}  {PARSER}\n")
    (DEST/"staging.json").write_text(json.dumps({
        "destination":str(DEST),"time":10,
        "expected_failed_checks":{"aligned":[],"jitter":["disk_ballistic_residual"]},
        "files_sha256":{name:digest(DEST/name) for name in sources},
        "external_parser_sha256":digest(PARSER),
        "scope":"Static geometry feasibility; failed native jitter gate remains failed"},indent=2)+"\n")
    print(json.dumps({"destination":str(DEST),"independent_cases":list(RECEIPTS)}))


if __name__=="__main__":
    main()

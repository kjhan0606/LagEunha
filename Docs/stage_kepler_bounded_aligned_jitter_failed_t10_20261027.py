"""Stage exploratory comparison of an accepted and an explicitly failed t10 state."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DEST = ROOT/"aligned_jitter_failed_t10_cpu_v2"
RUNNER = "kepler_bounded_aligned_jitter_failed_t10_cpu_20261027.slurm"
FILES = (
    RUNNER,
    "audit_kepler_bounded_radial_tails_20261027.py",
    "audit_kepler_bounded_aligned_jitter_pair_20261027.py",
    "plot_kepler_bounded_aligned_jitter_density_20261027.py",
    "plot_kepler_bounded_periodic_density_20261027.py",
    "plot_kepler_native_density_long_20261026.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_native_ballistic_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)
PARSER = Path("/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/audit_snapshot.py")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4*1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    aligned_path = ROOT/"to10_near32_a10_v1/audit_t10.json"
    jitter_path = ROOT/"jitter005_to10_a10_v1/audit_t10.json"
    receipts = []
    for path, expected in ((aligned_path,set()),
                           (jitter_path,{"disk_ballistic_residual"})):
        data=json.loads(path.read_text())
        failed={k for k,v in data["checks"].items() if not v}
        if (data["time"] != 10 or failed != expected
                or digest(Path(data["initial"])) != data["initial_sha256"]
                or digest(Path(data["endpoint"])) != data["endpoint_sha256"]):
            raise ValueError(f"Unexpected or changed native t10 result: {path}")
        receipts.append({"path":str(path),"sha256":digest(path),
                         "endpoint_sha256":data["endpoint_sha256"],
                         "failed_checks":sorted(failed)})
    subprocess.run(["bash","-n",str(DOCS/RUNNER)],check=True)
    if not PARSER.is_file():
        raise FileNotFoundError(PARSER)
    DEST.mkdir()
    for name in FILES:
        shutil.copy2(DOCS/name,DEST/name)
    hashes={name:digest(DEST/name) for name in FILES}
    with (DEST/"prestaged.sha256").open("x") as stream:
        for name,sha in sorted(hashes.items()):
            stream.write(f"{sha}  {name}\n")
        stream.write(f"{digest(PARSER)}  {PARSER}\n")
    manifest={"destination":str(DEST),"runner":RUNNER,
              "native_receipts":receipts,
              "prestaged_sha256":hashes,
              "external_snapshot_parser":{"path":str(PARSER),"sha256":digest(PARSER)},
              "interpretation":"Exploratory diagnosis; jitter disk_ballistic_residual remains failed"}
    (DEST/"prestage_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps({"destination":str(DEST),"native_receipts":receipts}))


if __name__ == "__main__":
    main()

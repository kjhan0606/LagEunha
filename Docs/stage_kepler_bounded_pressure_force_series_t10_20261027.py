"""Stage the CPU-only saved-checkpoint force-growth audit."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


DOCS=Path(__file__).resolve().parent
ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DEST=ROOT/"pressure_force_series_t10_cpu_v1"
FILES=("kepler_bounded_pressure_force_series_t10_cpu_20261027.slurm",
       "audit_kepler_bounded_pressure_force_t10_20261027.py",
       "summarize_kepler_bounded_pressure_force_series_20261027.py",
       "audit_kepler_sparse_area_t0825_20261027.py",
       "laguerre_nonlinear_volume_probe_t6_20261026.py",
       "prepare_kepler_sparse_area_restart_20261027.py")
EXTERNAL=(Path("/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/audit_snapshot.py"),
          Path("/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/varpress_w0_geom_cfl_long_20261026/audit_varpress_checkpoint_20261026.py"))


def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    subprocess.run(["bash","-n",str(DOCS/FILES[0])],check=True)
    receipts=[]
    for subdir,expected in (("to10_near32_a10_v1",set()),
                            ("jitter005_to10_a10_v1",{"disk_ballistic_residual"})):
        path=ROOT/subdir/"audit_t10.json"
        r=json.loads(path.read_text())
        failed={k for k,v in r["checks"].items() if not v}
        if r["time"]!=10 or failed!=expected or digest(Path(r["endpoint"]))!=r["endpoint_sha256"]:
            raise ValueError(f"Changed native t10 receipt: {path}")
        receipts.append({"path":str(path),"sha256":digest(path),"failed_checks":sorted(failed)})
    DEST.mkdir()
    for name in FILES:
        shutil.copy2(DOCS/name,DEST/name)
    with (DEST/"prestaged.sha256").open("x") as stream:
        for name in FILES:
            stream.write(f"{digest(DEST/name)}  {name}\n")
        for path in EXTERNAL:
            stream.write(f"{digest(path)}  {path}\n")
    (DEST/"prestage_manifest.json").write_text(json.dumps({
        "destination":str(DEST),"native_receipts":receipts,
        "source_sha256":{name:digest(DEST/name) for name in FILES},
        "external_sha256":{str(path):digest(path) for path in EXTERNAL},
        "scope":"CPU-only diagnostic of saved full native checkpoints; no calculation restart"},indent=2)+"\n")
    print(json.dumps({"destination":str(DEST),"native_receipts":receipts}))


if __name__=="__main__":
    main()

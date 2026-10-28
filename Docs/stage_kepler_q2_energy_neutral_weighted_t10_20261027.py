"""Stage independent native weighted t=10 continuations after one-step gates."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS=Path(__file__).resolve().parent
ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PILOT=ROOT/"gpu_energy_neutral_weighted_t6501_a10_v2"
RUNNER="kepler_q2_energy_neutral_weighted_t10_a10_20261027.slurm"
PARAMS_SHA="2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"
HELPERS=(
    "audit_kepler_q2_weighted_t10_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "laguerre_nonlinear_volume_probe_t6_20261026.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)
PARENTS={
    "aligned":"KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json",
    "jitter":"KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json",
}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("destination",type=Path)
    args=ap.parse_args()
    dest=args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    subprocess.run(["bash","-n",str(DOCS/RUNNER)],check=True)
    pilot_manifest=json.loads((PILOT/"stage_manifest.json").read_text())
    binary_sha=pilot_manifest["binary_sha256"]
    accepted={}
    for label,receipt_name in PARENTS.items():
        parent=json.loads((DOCS/receipt_name).read_text())
        pilot=PILOT/label
        audit=json.loads((pilot/"audit_t6501.json").read_text())
        endpoint=Path(audit["endpoint"])
        initial=Path(parent["initial"])
        if (not all(audit["checks"].values()) or abs(audit["endpoint_time"]-6.501)>1e-10
                or endpoint.parent!=pilot or digest(endpoint)!=audit["endpoint_sha256"]
                or abs(snapshot_time(endpoint,55351)-6.501)>1e-10
                or audit["restart_sha256"]!=parent["endpoint_sha256"]
                or digest(initial)!=parent["initial_sha256"]
                or abs(snapshot_time(initial,55351))>1e-12
                or digest(pilot/"eunha2.exe")!=binary_sha
                or digest(pilot/"params.dat")!=PARAMS_SHA):
            raise ValueError(f"{label}: pilot or accepted parent changed")
        accepted[label]=(pilot,audit,initial)
    dest.mkdir(parents=True)
    manifest={"run":str(dest),"target_time":10.,"binary_sha256":binary_sha,
              "pilot_run":str(PILOT),"cases":{},
              "scope":"Independent weighted native t=10 continuations from accepted t=6.501; not GIZMO t=600."}
    for label,(pilot,audit,initial) in accepted.items():
        case=dest/label
        case.mkdir()
        sources={
            "eunha2.exe":pilot/"eunha2.exe",
            "params.dat":pilot/"params.dat",
            "khout.000001.dat":Path(audit["endpoint"]),
            "initial.dat":initial,
            "disk_active_mask.bin":PILOT/"disk_active_mask.bin",
            "job.slurm":DOCS/RUNNER,
            **{name:DOCS/name for name in HELPERS},
        }
        checksums={}
        for name,source in sources.items():
            shutil.copy2(source,case/name)
            checksums[name]=digest(case/name)
        with (case/"inputs.sha256").open("x") as stream:
            for name,value in sorted(checksums.items()):
                stream.write(f"{value}  {name}\n")
        manifest["cases"][label]={"parent_audit_sha256":digest(pilot/"audit_t6501.json"),
                                  "pilot_endpoint_sha256":audit["endpoint_sha256"],
                                  "initial_sha256":digest(initial),"inputs":checksums}
    with (dest/"stage_manifest.json").open("x") as stream:
        json.dump(manifest,stream,indent=2)
        stream.write("\n")
    print(json.dumps({"run":str(dest),"binary_sha256":binary_sha,
                      "pilot_endpoint_sha256":{label:item[1]["endpoint_sha256"]
                                               for label,item in accepted.items()}},indent=2))


if __name__=="__main__":
    main()

"""Stage a hash-gated fresh weighted jitter pilot without touching prior runs."""

import hashlib
import json
from pathlib import Path
import shutil


ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DOCS=Path(__file__).resolve().parent
SOURCE=ROOT/"jitter005_t001_a10_v1"
BUILD=ROOT/"gpu_energy_neutral_native_build_v3_jacobi_choice"
MASK=ROOT/"gpu_energy_neutral_weighted_to10_a10_v1"/"disk_active_mask.bin"
RUN=ROOT/"gpu_energy_neutral_weighted_fresh_jitter64_t001_a10_v1"
EXPECTED_INITIAL="8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b"
EXPECTED_BINARY="7358c54996b820a7b45bd53bcf6f788f125754b3158845d2b4265afead570f38"
EXPECTED_MASK="b560aca44c81e4560f8371c6d979e1afbe853f8cd18b592d0fc19a58becc8f4c"


def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def main():
    origins={
        "eunha2.exe":BUILD/"eunha2_energy_neutral_v3.exe",
        "glass.bin":SOURCE/"glass.bin",
        "areas.bin":SOURCE/"areas.bin",
        "params.dat":SOURCE/"params.dat",
        "accepted_initial.dat":SOURCE/"khout.000000.dat",
        "ic_manifest.json":SOURCE/"ic_manifest.json",
        "disk_active_mask.bin":MASK,
        "job.slurm":DOCS/"kepler_q2_weighted_fresh_jitter64_t001_20261027.slurm",
    }
    for name in (
        "audit_kepler_q2_weighted_fresh_pilot_20261027.py",
        "audit_kepler_sparse_area_t0825_20261027.py",
        "audit_kepler_sparse_state_consistency_20261027.py",
        "laguerre_nonlinear_volume_probe_t6_20261026.py",
        "prepare_kepler_sparse_area_restart_20261027.py",
    ):
        origins[name]=DOCS/name
    expected={"accepted_initial.dat":EXPECTED_INITIAL,
              "eunha2.exe":EXPECTED_BINARY,"disk_active_mask.bin":EXPECTED_MASK}
    for name,value in expected.items():
        actual=digest(origins[name])
        if actual!=value:
            raise ValueError(f"{name} hash changed: {actual}")
    if RUN.exists():
        raise FileExistsError(RUN)
    RUN.mkdir()
    hashes={}
    for name,source in origins.items():
        target=RUN/name
        shutil.copy2(source,target)
        hashes[name]=digest(target)
        if hashes[name]!=digest(source):
            raise ValueError(f"Copy changed: {name}")
    with (RUN/"inputs.sha256").open("x") as output:
        for name,value in sorted(hashes.items()):
            output.write(f"{value}  {name}\n")
    manifest={
        "run":str(RUN),"source_preflight":str(SOURCE),"weighted_binary":str(BUILD),
        "accepted_initial_sha256":EXPECTED_INITIAL,"binary_sha256":EXPECTED_BINARY,
        "mask_sha256":EXPECTED_MASK,"target_time":.001,"jacobi_iterations":64,
        "inputs":{name:{"source":str(origins[name]),"sha256":hashes[name]}
                  for name in sorted(hashes)},
        "scope":"Fresh t=0 experimental weighted jitter pilot; prerequisite for a "
                "new t=0→6.5 chain, not a restart of the prior W0 history.",
    }
    with (RUN/"stage_manifest.json").open("x") as output:
        json.dump(manifest,output,indent=2)
        output.write("\n")
    print(json.dumps({"run":str(RUN),"binary_sha256":EXPECTED_BINARY,
                      "input_count":len(hashes)},indent=2))


if __name__=="__main__":
    main()

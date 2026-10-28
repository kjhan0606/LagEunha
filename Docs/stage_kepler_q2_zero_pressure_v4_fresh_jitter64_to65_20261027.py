"""Stage the t=.001→6.5 continuation from verified v4 pilot 437452."""

import hashlib
import json
from pathlib import Path
import shutil


ROOT=Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DOCS=Path(__file__).resolve().parent
PILOT=ROOT/"gpu_energy_neutral_zero_pressure_v4_t001_a10_v1"
RUN=ROOT/"gpu_energy_neutral_zero_pressure_v4_fresh_jitter64_to65_a10_v1"
EXPECTED_INITIAL="8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b"
EXPECTED_BINARY="3120aea6dad944b144e04d4597cb3f644abc6c803ff2318b5862b6bf8eec0adc"
EXPECTED_MASK="b560aca44c81e4560f8371c6d979e1afbe853f8cd18b592d0fc19a58becc8f4c"
ZERO_GPU=ROOT/'zero_pressure_gradient_gpu_test_v1/fixed_result.log'
ZERO_GPU_SHA='b32e0de19ab0b9b3d368f54b9b1fe4afd99cdfc76467cc750e122a1b27601a9a'


def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def main():
    pilot=json.load(open(PILOT/'audit_t001.json'))
    pair=json.load(open(PILOT/'pair_audit.json'))
    if (not all(pilot['checks'].values()) or not all(pair['checks'].values())
            or pilot['binary_sha256']!=EXPECTED_BINARY
            or pilot['reference_sha256']!=EXPECTED_INITIAL
            or pair['v4_endpoint_sha256']!=pilot['endpoint_sha256']
            or digest(Path(pilot['endpoint']))!=pilot['endpoint_sha256']):
        raise ValueError('The v4 pilot, pair gate or endpoint did not pass')
    if digest(ZERO_GPU)!=ZERO_GPU_SHA:
        raise ValueError('The exact-zero CUDA regression result changed')
    origins={
        "eunha2.exe":PILOT/"eunha2.exe",
        "params.dat":PILOT/"params.dat",
        "initial.dat":PILOT/"accepted_initial.dat",
        "disk_active_mask.bin":PILOT/"disk_active_mask.bin",
        "job.slurm":DOCS/"kepler_q2_zero_pressure_v4_fresh_jitter64_to65_20261027.slurm",
    }
    for name in (
        "audit_kepler_q2_weighted_equilibrium_20261027.py",
        "audit_kepler_bounded_native_evolution_20261027.py",
        "audit_kepler_sparse_native_ballistic_20261027.py",
        "audit_kepler_sparse_area_t0825_20261027.py",
        "audit_kepler_sparse_state_consistency_20261027.py",
        "laguerre_nonlinear_volume_probe_t6_20261026.py",
        "prepare_kepler_sparse_area_restart_20261027.py",
    ):
        origins[name]=DOCS/name
    expected={"initial.dat":EXPECTED_INITIAL,
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
        "run":str(RUN),"pilot":str(PILOT),"dependency":"afterok:437452",
        "accepted_initial_sha256":EXPECTED_INITIAL,"binary_sha256":EXPECTED_BINARY,
        "zero_gradient_gpu_result":str(ZERO_GPU),"zero_gradient_gpu_sha256":ZERO_GPU_SHA,
        "mask_sha256":EXPECTED_MASK,"start_time":.001,"target_time":6.5,
        "jacobi_iterations":64,
        "inputs":{name:{"source":str(origins[name]),"sha256":hashes[name]}
                  for name in sorted(hashes)},
        "scope":"Conditional fresh weighted jitter continuation. Runtime accepts only "
                "the hash-verified complete passing pilot endpoint.",
    }
    with (RUN/"stage_manifest.json").open("x") as output:
        json.dump(manifest,output,indent=2)
        output.write("\n")
    print(json.dumps({"run":str(RUN),"dependency":manifest["dependency"],
                      "input_count":len(hashes)},indent=2))


if __name__=="__main__":
    main()

#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
source "$root/run_env.sh"
export OMP_NUM_THREADS=1 I_MPI_HYDRA_BOOTSTRAP=fork NATIVE_MPI_RANKS=1
export EUNHA_IC=affine_probe HYDRO_TSTOP=0.000000000001 EUNHA_DUMP_DT=0.000000001 GFS_GPU_FACE_DUMP_CALL=1 KEPLER_RESTART_STEP=1 GFS_GPU_FINAL_VOLUME=1
unset I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
python3 - <<'PY'
from pathlib import Path
import hashlib,json
root=Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
base=root/'irregular_affine_probe'
binary=root/'affine_limiter_intervention/eunha2_affine_probe.exe'
manifest=json.loads((base/'manifest.json').read_text())
cases=[('rigid_irregular_lsf_default','rigid',0,0),
       ('rigid_irregular_lsf_unlimited','rigid',1,0),
       ('rigid_irregular_lsf_unmerged','rigid',0,1),
       ('rigid_irregular_lsf_unlimited_unmerged','rigid',1,1),
       ('strain_irregular_lsf_unlimited_unmerged','strain',1,1)]
for name,flow,v,m in cases:
    dest=base/name
    dest.mkdir(exist_ok=False)
    original=base/f'{flow}_irregular_lsf'
    (dest/'khout.000001.dat').symlink_to(original/'khout.000001.dat')
    (dest/'params.dat').write_text((original/'params.dat').read_text())
    (dest/'eunha2.exe').symlink_to(binary)
    (dest/'OWNER.txt').write_text('Codex affine causal diagnostic; native source controls restricted to affine_probe\n')
    manifest['cases'].append({'name':name,'flow':flow,'gradient_method':1,'OoA':.333333,
                             'unlimited_velocity':v,'unmerged_faces':m,
                             'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest()})
(base/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
PY
for variant in rigid_irregular_lsf_default rigid_irregular_lsf_unlimited rigid_irregular_lsf_unmerged rigid_irregular_lsf_unlimited_unmerged strain_irregular_lsf_unlimited_unmerged; do
 export GFS_AFFINE_UNLIMITED_VELOCITY=0 GFS_AFFINE_UNMERGED_FACES=0
 case "$variant" in *unlimited*) export GFS_AFFINE_UNLIMITED_VELOCITY=1;; esac
 case "$variant" in *unmerged*) export GFS_AFFINE_UNMERGED_FACES=1;; esac
 cd "$root/irregular_affine_probe/$variant"
 test ! -e log
 python3 - <<'PY'
from pathlib import Path
import os
keys=('OMP_NUM_THREADS','NATIVE_MPI_RANKS','EUNHA_IC','HYDRO_TSTOP','EUNHA_DUMP_DT',
      'GFS_GPU_FACE_DUMP_CALL','KEPLER_RESTART_STEP','GFS_GPU_FINAL_VOLUME',
      'GFS_AFFINE_UNLIMITED_VELOCITY','GFS_AFFINE_UNMERGED_FACES')
Path('environment.txt').write_text(''.join(f'{k}={os.environ.get(k, "")}\n' for k in keys))
PY
 sha256sum ./eunha2.exe params.dat khout.000001.dat > inputs.sha256
 python3 "$root/run_native.py"
done
python3 "$root/analyze_irregular_affine.py"

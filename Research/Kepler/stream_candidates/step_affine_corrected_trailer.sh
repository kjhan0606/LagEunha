#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
source "$root/run_env.sh"
export OMP_NUM_THREADS=1 I_MPI_HYDRA_BOOTSTRAP=fork NATIVE_MPI_RANKS=1
export EUNHA_IC=affine_probe HYDRO_TSTOP=0.000000000001 EUNHA_DUMP_DT=0.000000001 GFS_GPU_FACE_DUMP_CALL=1 KEPLER_RESTART_STEP=1 GFS_GPU_FINAL_VOLUME=1
unset I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
python3 - <<'PY'
from pathlib import Path
import hashlib,json,struct
root=Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
base=root/'irregular_affine_probe_v2'
base.mkdir(exist_ok=False)
binary=root/'affine_limiter_intervention/eunha2_affine_probe.exe'
old=root/'irregular_affine_probe'
manifest=json.loads((old/'manifest.json').read_text())
manifest['cases']=[]
manifest['correction']='v1 custom trailer accidentally wrote gamma into GAS_dtold. Particle payload is unchanged here; trailer corrected to standard cold start t0, dt1e-7, dtold1e-7. Old outputs preserved. Initial gradient/geometry findings are independent of this field, pair pressure caps are not.'
cases=[('rigid_limited_merged','rigid',0,0),
       ('rigid_unlimited_merged','rigid',1,0),
       ('rigid_limited_unmerged','rigid',0,1),
       ('rigid_unlimited_unmerged','rigid',1,1),
       ('strain_unlimited_unmerged','strain',1,1)]
for name,flow,v,m in cases:
    dest=base/name
    dest.mkdir()
    original=old/f'{flow}_irregular_lsf'
    src=(original/'khout.000001.dat').read_bytes()
    payload=src[:-52]
    trailer=struct.pack('<ddfiiddff',0.,1e-7,1e-7,256,256,24.,24.,0.,0.)
    (dest/'khout.000001.dat').write_bytes(payload+trailer)
    (dest/'params.dat').write_text((original/'params.dat').read_text())
    (dest/'eunha2.exe').symlink_to(binary)
    (dest/'OWNER.txt').write_text('Codex affine causal diagnostic; correct native KH trailer; one native step\n')
    manifest['cases'].append({'name':name,'flow':flow,'gradient_method':1,'OoA':.333333,
                             'unlimited_velocity':v,'unmerged_faces':m,
                             'initial_payload_sha256':hashlib.sha256(payload).hexdigest(),
                             'initial_sha256':hashlib.sha256(payload+trailer).hexdigest(),
                             'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest()})
(base/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
PY
for variant in rigid_limited_merged rigid_unlimited_merged rigid_limited_unmerged rigid_unlimited_unmerged strain_unlimited_unmerged; do
 export GFS_AFFINE_UNLIMITED_VELOCITY=0 GFS_AFFINE_UNMERGED_FACES=0
 case "$variant" in *unlimited*) export GFS_AFFINE_UNLIMITED_VELOCITY=1;; esac
 case "$variant" in *unmerged*) export GFS_AFFINE_UNMERGED_FACES=1;; esac
 cd "$root/irregular_affine_probe_v2/$variant"
 test ! -e log
 sha256sum ./eunha2.exe params.dat khout.000001.dat > inputs.sha256
 python3 "$root/run_native.py"
done
python3 "$root/analyze_irregular_affine.py" irregular_affine_probe_v2

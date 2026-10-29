"""Replay an existing Q2 pre-impact checkpoint with experimental thermal64."""
import hashlib
import json
from pathlib import Path
import shutil

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
parent=Path('/gpfs/kjhan/LagForce/Kepler_Q2_WallCrossing_Replay_20261029/native_events_v1/crossing_to484')
root=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallImpact_20261029')
src=(parent/'khout.011580.dat').resolve()
exe=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallCandidate_20261029/project/eunha2.exe')
if sha(src)!='c4b11f1996349f5ce8ab664569fa4415cd91b5f836fe9579c98f7dcd51aff758':raise ValueError('Unexpected preimpact state')
if sha(exe)!='ff3ae9b86c65e86a28dfe246a33186e9387361ee54c83439bb17bc52f85e317a':raise ValueError('Unexpected candidate')
if root.exists():raise FileExistsError(root)
root.mkdir();(root/src.name).symlink_to(src);(root/'eunha2.exe').symlink_to(exe)
shutil.copy2(parent/'params.dat',root/'params.dat')
run=(parent/'run.sh').read_text().replace(str(parent),str(root)).replace('HYDRO_TSTOP=48.4','HYDRO_TSTOP=48.378')
old='export GFS_RK4_AUDIT_STEP_MIN=0 GFS_RK4_AUDIT_STEP_MAX=1000000\nexport GFS_RK4_AUDIT_STAGE=k1_rhs GFS_RK4_AUDIT_IDS=55283'
if run.count(old)!=1:raise ValueError('Unexpected diagnostics template')
run=run.replace(old,'unset GFS_RK4_AUDIT_STEP_MIN GFS_RK4_AUDIT_STEP_MAX GFS_RK4_AUDIT_STAGE GFS_RK4_AUDIT_IDS')
run=run.replace('export GFS_KEPLER_WALL_EVENTS=1 GFS_GPU_TESS_DUMP=0',
    'export GFS_KEPLER_WALL_EVENTS=1 GFS_GPU_TESS_DUMP=0 GFS_THERMAL64=1 GFS_NATIVE_CONST_CONFIG=0\nunset GFS_THERMAL64_INIT_ONLY GFS_RK4_DIRECT_U GFS_GPU_FACE_DUMP_FIRST GFS_GPU_FACE_DUMP_LAST\nulimit -f 131072')
(root/'run.sh').write_text(run)
stage=json.loads((parent/'stage.json').read_text());stage.update(target_time=48.378,binary_sha256=sha(exe),
    role='Experimental thermal64 collision replay from accepted Q2 checkpoint; no Q1 stability claim',
    output_policy='One final snapshot (~32 MB), no particle TSV')
(root/'stage.json').write_text(json.dumps(stage,indent=2)+'\n')
(root/'job_ada.sh').write_text('''#!/bin/bash
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 CUDA_VISIBLE_DEVICES=GPU-42b3367f-1df7-2b2d-7961-7781c87b022d
bash /gpfs/kjhan/LagForce/Kepler_Thermal64WallImpact_20261029/run.sh
''')
print(json.dumps(stage,indent=2))

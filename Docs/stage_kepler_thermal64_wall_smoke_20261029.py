"""Stage a very short forced-rollback smoke from an accepted bounded checkpoint."""
import hashlib
from pathlib import Path
import shutil

root = Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallSmoke_20261029')
parent = Path('/gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029/to0p4125')
initial = Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/gravity_cap_preflight_v1/native_eta006_ada/khout.000022.dat')
exe = Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallCandidate_20261029/project/eunha2.exe')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()
if sha(initial) != '0f60d80dec2668a7e53a6e99ec84760af824ed70fb7a868ddd605fe53ffe11bc': raise ValueError('Unexpected accepted start')
if sha(exe) != 'ff3ae9b86c65e86a28dfe246a33186e9387361ee54c83439bb17bc52f85e317a': raise ValueError('Unexpected candidate binary')
if root.exists(): raise FileExistsError(root)
root.mkdir()
(root/'khout.000022.dat').symlink_to(initial)
(root/'eunha2.exe').symlink_to(exe)
shutil.copy2(parent/'params.dat', root/'params.dat')
text = (parent/'run.sh').read_text()
text = text.replace('cd '+str(parent), 'cd '+str(root))
text = text.replace('OMP_NUM_THREADS=8', 'OMP_NUM_THREADS=2')
text = text.replace('HYDRO_TSTOP=0.41249999999999998', 'HYDRO_TSTOP=0.0016')
line = '/home/kjhan/miniconda3/bin/python3 /gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/post_t10_restart_dt_h200_v1/run_kepler_q1_restart_dt_bundle_20261029.py validate-environment /gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029/to0p4125/case.json'
if text.count(line) != 1: raise ValueError('Unexpected template')
text = text.replace(line, 'export GFS_THERMAL64=1 GFS_WALL_FORCE_REJECT_ONCE=1\nunset GFS_THERMAL64_INIT_ONLY GFS_RK4_DIRECT_U\nulimit -f 131072')
(root/'run.sh').write_text(text)
(root/'job_ada.sh').write_text('''#!/bin/bash
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 CUDA_VISIBLE_DEVICES=GPU-42b3367f-1df7-2b2d-7961-7781c87b022d
bash /gpfs/kjhan/LagForce/Kepler_Thermal64WallSmoke_20261029/run.sh 22
''')
print('staged', root, 'from accepted t=0.0015284035789018555 to t=0.0016 with one forced rejection')

"""Stage an early fixed-mass W0 thermal64 Q1 run from the accepted eta=.06 restart."""
import hashlib
import json
from pathlib import Path
import shutil

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Early_H200_20261029')
base=Path('/gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029/to0p4125')
seed=(base/'khout.000022.dat').resolve()
exe=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallCandidate_20261029/project/eunha2.exe')
if sha(seed)!='0f60d80dec2668a7e53a6e99ec84760af824ed70fb7a868ddd605fe53ffe11bc':raise ValueError('Changed accepted restart')
if sha(exe)!='ff3ae9b86c65e86a28dfe246a33186e9387361ee54c83439bb17bc52f85e317a':raise ValueError('Changed thermal candidate')
screen=json.loads((base/'core_screen.json').read_text());gate=json.loads((base/'checkpoint_gate.json').read_text())
if not screen['passed'] or len(gate['checks'])!=28 or not all(gate['checks'].values()):raise ValueError('Missing certified corresponding legacy checkpoint')
if root.exists():raise FileExistsError(root)
root.mkdir();(root/seed.name).symlink_to(seed);(root/'eunha2.exe').symlink_to(exe);shutil.copy2(base/'params.dat',root/'params.dat')
run=(base/'run.sh').read_text().replace(str(base),str(root)).replace('unset I_MPI_HYDRA_BOOTSTRAP','export I_MPI_HYDRA_BOOTSTRAP=slurm')
validator='/home/kjhan/miniconda3/bin/python3 /gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/post_t10_restart_dt_h200_v1/run_kepler_q1_restart_dt_bundle_20261029.py validate-environment '+str(root)+'/case.json'
if run.count(validator)!=1:raise ValueError('Unexpected original runner')
run=run.replace(validator,'export GFS_THERMAL64=1\nunset GFS_THERMAL64_INIT_ONLY GFS_RK4_DIRECT_U\nulimit -f 131072\nsha256sum -c inputs.sha256')
(root/'run.sh').write_text(run)
files=[seed,exe,root/'params.dat',root/'run.sh']
(root/'inputs.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in files))
info={'source':str(seed),'source_sha256':sha(seed),'start_time':0.0015284035789018555,'target_time':0.4125,
      'eta':0.06,'binary':str(exe),'binary_sha256':sha(exe),'legacy_reference':str(base/'khout.005400.dat'),
      'legacy_reference_sha256':sha(base/'khout.005400.dat'),'reference_core_screen_sha256':sha(base/'core_screen.json'),
      'output_policy':'One 52 MB final snapshot, no intermediate dumps or full-particle TSV.',
      'scope':'Early native thermal64 Q1 accuracy comparison, not a fresh t0 run or t6.5/long-time proof.'}
(root/'contract.json').write_text(json.dumps(info,indent=2)+'\n')
job='''#!/bin/bash
#SBATCH --job-name=KepQ1ThermEarly
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --chdir='''+str(root)+'''
#SBATCH --output='''+str(root)+'''/slurm-%j.out
#SBATCH --error='''+str(root)+'''/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=8 OMP_STACKSIZE=256M OPENBLAS_NUM_THREADS=1
hostname
nvidia-smi --query-gpu=name,uuid --format=csv,noheader
bash run.sh 22
'''
(root/'job.slurm').write_text(job)
print(json.dumps(info,indent=2))

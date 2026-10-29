"""Stage same-parent t10 -> t10.1 thermal64 comparison with the existing half control."""
import hashlib
import json
from pathlib import Path
import shutil

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

base=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/post_t10_restart_dt_h200_v1/controls/half')
root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Half_H200_20261029')
state=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/chain_h200_v2/to10p0/khout.020790.dat')
binary=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallCandidate_20261029/project/eunha2.exe')
if sha(state)!='58afbe3ab219ddef0b8dbb9b032789236ec01815bcf80457a5d73e6b90bd6c28':raise ValueError('Unexpected t10 checkpoint')
if sha(binary)!='ff3ae9b86c65e86a28dfe246a33186e9387361ee54c83439bb17bc52f85e317a':raise ValueError('Unexpected candidate')
case=json.loads((base/'case.json').read_text())
if case['target']!=10.1 or case['source']['sha256']!=sha(state) or case['cap']!=4.760027150041424e-05:raise ValueError('Unexpected single-precision comparison contract')
if root.exists():raise FileExistsError(root)
root.mkdir();(root/state.name).symlink_to(state);(root/'eunha2.exe').symlink_to(binary);shutil.copy2(base/'params.dat',root/'params.dat')
run=(base/'run.sh').read_text().replace(str(base),str(root))
run=run.replace('export GFS_RK4_GAP_PROFILE=1 GFS_NATIVE_CONST_CONFIG=1','export GFS_RK4_GAP_PROFILE=1 GFS_NATIVE_CONST_CONFIG=0')
old='/home/kjhan/miniconda3/bin/python3 /gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/post_t10_restart_dt_h200_v1/run_kepler_q1_restart_dt_bundle_20261029.py validate-environment '+str(root)+'/case.json'
if run.count(old)!=1:raise ValueError('Unexpected control runner')
run=run.replace(old,'export GFS_THERMAL64=1\nunset GFS_THERMAL64_INIT_ONLY GFS_RK4_DIRECT_U\nulimit -f 131072\nsha256sum -c inputs.sha256')
(root/'run.sh').write_text(run)
files=[state,binary,root/'params.dat',root/'run.sh']
(root/'inputs.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in files))
contract={'source':str(state),'source_sha256':sha(state),'binary':str(binary),'binary_sha256':sha(binary),
          'reference_half':str(base),'target':10.1,'dt_cap':case['cap'],
          'output_policy':'One final 52 MB snapshot; no intermediate/full-particle TSV; preserve all inputs.',
          'scope':'Compare arithmetic over t10 to 10.1 from the same already-damaged state; not a t>10 stability or core-repair claim.'}
(root/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
job='''#!/bin/bash
#SBATCH --job-name=KepQ1Therm64Half
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --chdir='''+str(root)+'''
#SBATCH --output='''+str(root)+'''/slurm-%j.out
#SBATCH --error='''+str(root)+'''/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=8 OMP_STACKSIZE=256M OPENBLAS_NUM_THREADS=1
hostname
nvidia-smi --query-gpu=name,uuid --format=csv,noheader
bash run.sh 20790
'''
(root/'job.slurm').write_text(job)
print(json.dumps(contract,indent=2))

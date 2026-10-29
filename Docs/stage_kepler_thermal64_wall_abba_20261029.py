"""Stage same-state ABBA wall-hoist parity and timing on one H200 allocation."""
import hashlib
import json
from pathlib import Path
import shutil

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
root=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallABBA_H200_20261029')
parent=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Half_H200_20261029')
seed=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/chain_h200_v2/to10p0/khout.020790.dat')
a=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallCandidate_20261029/project/eunha2.exe')
b=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallHoist_20261029/project/eunha2.exe')
expected={'seed':'58afbe3ab219ddef0b8dbb9b032789236ec01815bcf80457a5d73e6b90bd6c28','A':'ff3ae9b86c65e86a28dfe246a33186e9387361ee54c83439bb17bc52f85e317a','B':'9ae8ddd03bf48cda95352ef9636ec231ae9efc4433c44f8015b5d020adeb3610'}
if any(sha(p)!=expected[k] for k,p in [('seed',seed),('A',a),('B',b)]):raise ValueError('Changed frozen source or binary')
if root.exists():raise FileExistsError(root)
root.mkdir();entries=[]
for name,exe in [('A1',a),('B1',b),('B2',b),('A2',a)]:
    d=root/name;d.mkdir();(d/seed.name).symlink_to(seed);(d/'eunha2.exe').symlink_to(exe)
    shutil.copy2(parent/'params.dat',d/'params.dat')
    script=(parent/'run.sh').read_text().replace(str(parent),str(d)).replace('HYDRO_TSTOP=10.1','HYDRO_TSTOP=10.001')
    if script.count('HYDRO_TSTOP=10.001')!=1:raise ValueError('Unexpected native target')
    (d/'run.sh').write_text(script)
    files=[seed,exe,d/'params.dat',d/'run.sh']
    (d/'inputs.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in files))
    entries.append({'case':name,'binary_sha256':sha(exe),'run_sha256':sha(d/'run.sh'),'inputs_sha256':sha(d/'inputs.sha256')})
contract={'seed_sha256':sha(seed),'target_time':10.001,'step_cap':4.760027150041424e-05,'order':entries,
          'outputs':'Four ~52 MB native endpoints and small logs, no TSV; no overwrite.',
          'scope':'One-node same-parent short ABBA for physics parity and step timing; no long-time stability or MPI scaling claim.'}
(root/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
job='''#!/bin/bash
#SBATCH --job-name=KepThermWallABBA
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --chdir='''+str(root)+'''
#SBATCH --output='''+str(root)+'''/slurm-%j.out
#SBATCH --error='''+str(root)+'''/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=8 OMP_STACKSIZE=256M OPENBLAS_NUM_THREADS=1
hostname
nvidia-smi --query-gpu=name,uuid --format=csv,noheader
for case in A1 B1 B2 A2; do
  cd '''+str(root)+'''/"$case"
  bash run.sh 20790
done
'''
(root/'job.slurm').write_text(job)
print(json.dumps(contract,indent=2))

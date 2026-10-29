"""Same-H200 zero-time thermal initialization and uncapped t10 control."""
import hashlib
import json
from pathlib import Path
import shutil

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Coarse_H200_20261029')
half=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Half_H200_20261029')
contract=json.loads((half/'contract.json').read_text());seed=Path(contract['source']);exe=Path(contract['binary'])
if sha(seed)!=contract['source_sha256'] or sha(exe)!=contract['binary_sha256']:raise ValueError('Changed source or binary')
comparison=json.loads((half/'arithmetic_comparison.json').read_text())
if not comparison['passed'] or comparison['source_sha256']!=sha(seed):raise ValueError('Half control not verified')
if root.exists():raise FileExistsError(root)
root.mkdir();cases=[]
for name in ('init','coarse'):
    d=root/name;d.mkdir();(d/seed.name).symlink_to(seed);(d/'eunha2.exe').symlink_to(exe)
    shutil.copy2(half/'params.dat',d/'params.dat')
    run=(half/'run.sh').read_text().replace(str(half),str(d))
    cap='export GFS_DIAG_DT_MAX=4.7600271500414237e-05\n'
    if run.count(cap)!=1:raise ValueError('Unexpected cap in half runner')
    run=run.replace(cap,'')
    if name=='init':
        marker='export GFS_THERMAL64=1\n'
        if run.count(marker)!=1:raise ValueError('Unexpected thermal config')
        run=run.replace(marker,marker+'export GFS_THERMAL64_INIT_ONLY=1 KEPLER_RESTART_STEP=20790 GFS_THERMAL64_INIT_TIME=10\n')
        run=run.replace('unset GFS_THERMAL64_INIT_ONLY GFS_RK4_DIRECT_U','unset GFS_RK4_DIRECT_U')
        if 'unset GFS_RK4_FAILURE_SNAPSHOT GFS_DIAG_CHECKPOINT_STEP' not in run:raise ValueError('Unexpected snapshot env')
        run=run.replace('unset GFS_RK4_FAILURE_SNAPSHOT GFS_DIAG_CHECKPOINT_STEP','unset GFS_DIAG_CHECKPOINT_STEP\nexport GFS_RK4_FAILURE_SNAPSHOT=1')
    (d/'run.sh').write_text(run)
    files=[seed,exe,d/'params.dat',d/'run.sh']
    (d/'inputs.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in files))
    cases.append({'name':name,'run_sha256':sha(d/'run.sh'),'inputs_sha256':sha(d/'inputs.sha256')})
info={'source_sha256':sha(seed),'binary_sha256':sha(exe),'cases':cases,'coarse_target':10.1,
      'output_policy':'One 52 MB zero-time diagnostic and one 52 MB final checkpoint, no full-particle TSV.',
      'scope':'Same-H200 initial volume reference then independent uncapped t10-to-10.1 evolution, compared with existing same-parent controls. Not long-time proof.'}
(root/'contract.json').write_text(json.dumps(info,indent=2)+'\n')
job='''#!/bin/bash
#SBATCH --job-name=KepQ1ThermCoarse
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
cd '''+str(root)+'''/init
bash run.sh 20790
cd '''+str(root)+'''/coarse
bash run.sh 20790
'''
(root/'job.slurm').write_text(job)
print(json.dumps(info,indent=2))

"""Stage eta=.06 thermal64 Ada chain from an accepted native checkpoint."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

docs=Path(__file__).resolve().parent
original=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Early_H200_20261029')
parent=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Early_Ada_20261029')
child=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64EarlyTo0825_Ada_20261029')
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def main():
    if parent.exists() or child.exists():raise FileExistsError('Ada target exists')
    cfg=json.loads((original/'contract.json').read_text())
    for p,k in [(Path(cfg['source']),'source_sha256'),(Path(cfg['binary']),'binary_sha256'),(Path(cfg['legacy_reference']),'legacy_reference_sha256')]:
        if sha(p)!=cfg[k]:raise ValueError('Changed input: '+str(p))
    if cfg['eta']!=.06 or cfg['target_time']!=.4125:raise ValueError('Unexpected original contract')
    parent.mkdir();(parent/Path(cfg['source']).name).symlink_to(Path(cfg['source']));(parent/'eunha2.exe').symlink_to(Path(cfg['binary']))
    shutil.copy2(original/'params.dat',parent/'params.dat')
    cfg['execution_host']='lageunha';cfg['cancelled_duplicate_h200_job']=439903
    cfg['scope']='Native thermal64 eta=.06 on Ada, from saved early checkpoint; not long-time proof.'
    (parent/'contract.json').write_text(json.dumps(cfg,indent=2)+'\n')
    script=(original/'run.sh').read_text()
    if script.count(str(original))!=1 or script.count('export I_MPI_HYDRA_BOOTSTRAP=slurm')!=1:raise ValueError('Unexpected source runner')
    script=script.replace(str(original),str(parent)).replace('export I_MPI_HYDRA_BOOTSTRAP=slurm','unset I_MPI_HYDRA_BOOTSTRAP')
    (parent/'run.sh').write_text(script)
    subprocess.run(['bash','-n',str(parent/'run.sh')],check=True)
    inputs=[Path(cfg['source']),Path(cfg['binary']),Path(cfg['legacy_reference']),parent/'params.dat',parent/'contract.json',parent/'run.sh']
    (parent/'inputs.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in inputs))
    subprocess.run(['/home/kjhan/miniconda3/bin/python3',str(docs/'run_kepler_q1_thermal64_ada_continuation_20261029.py'),'stage'],check=True)
    wrapper='''#!/bin/bash
set -eo pipefail
source /etc/profile.d/lmod.sh
set -u
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=8 OMP_STACKSIZE=256M OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH='''+str(docs)+'''
test "$(hostname | cut -d. -f1 | tr '[:upper:]' '[:lower:]')" = lageunha
nvidia-smi --query-gpu=name,uuid --format=csv,noheader
bash '''+str(parent)+'''/run.sh 22
cd '''+str(child)+'''
/home/kjhan/miniconda3/bin/python3 driver.py run
'''
    (parent/'chain.sh').write_text(wrapper);subprocess.run(['bash','-n',str(parent/'chain.sh')],check=True)
    print('Ada parent and gated continuation staged:',parent,child)
if __name__=='__main__':main()

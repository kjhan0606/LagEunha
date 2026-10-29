"""Stage one independent H200 ABBA benchmark; refuses existing output."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

docs=Path(__file__).resolve().parent
root=Path('/gpfs/kjhan/LagForce/Kepler_WallGuardOMP_ABBA_H200_20261029')
candidate=Path('/gpfs/kjhan/LagForce/Kepler_WallGuardOMP_20261029')
initial=(candidate/'omp/khout.000022.dat').resolve();binary=(candidate/'project/eunha2.exe').resolve()
params=candidate/'omp/params.dat';template=(candidate/'omp/run.sh').read_text()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if root.exists():raise FileExistsError(root)
if sha(initial)!='0f60d80dec2668a7e53a6e99ec84760af824ed70fb7a868ddd605fe53ffe11bc' or sha(binary)!='da70844984a52fec77e2d6dec451e9b5695e761ba93a42d0f0d079c1091edf44':raise ValueError('changed frozen input')
root.mkdir()
order=['A1','B1','B2','A2']
stage={'initial':str(initial),'initial_sha256':sha(initial),'binary_sha256':sha(binary),
       'target_time':0.003056807157803711,'order':order,'minimum_speedup':1.05,
       'output_policy':'4 final snapshots, ~208 MB total, no intermediate or full-particle TSV',
       'scope':'One allocated H200, one rank, eight threads, same binary and thermal64 initial state; guard off/on ABBA; no other jobs are dependencies.'}
(root/'stage.json').write_text(json.dumps(stage,indent=2)+'\n')
paths=[initial,binary,params,root/'stage.json']
for name in order:
    d=root/name;d.mkdir();(d/initial.name).symlink_to(initial);(d/'eunha2.exe').symlink_to(binary);shutil.copyfile(params,d/'params.dat')
    flag='1' if name.startswith('B') else '0'
    run=template.replace(str(candidate/'omp'),str(d)).replace('HYDRO_TSTOP=0.0018340842946822263','HYDRO_TSTOP=0.003056807157803711').replace('GFS_WALL_OMP_VALIDATE=1','GFS_WALL_OMP_VALIDATE='+flag)
    if str(candidate/'omp') in run or 'HYDRO_TSTOP=0.003056807157803711' not in run:raise ValueError('runner rewrite failed')
    (d/'run.sh').write_text(run);subprocess.run(['bash','-n',str(d/'run.sh')],check=True)
    files=[initial,binary,d/'params.dat',d/'run.sh'];(d/'inputs.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in files))
    paths.extend((d/'params.dat',d/'run.sh',d/'inputs.sha256'))
audits=[docs/x for x in ('audit_wall_guard_omp_abba_20261029.py','audit_kepler_restart_bytes_20261029.py',
        'audit_kepler_wall_events_interleaved_20261029.py','audit_kepler_native_radial_20261026.py','kepler_native_global_ids_20261029.py')]
paths.extend(audits);(root/'inputs.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in paths))
job=f'''#!/bin/bash
#SBATCH --job-name=KepWallOMP
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --chdir={root}
#SBATCH --output={root}/slurm-%j.out
#SBATCH --error={root}/slurm-%j.err
set -euo pipefail
source /etc/profile.d/lmod.sh
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export SLURM_EXPORT_ENV=ALL OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 OMP_PROC_BIND=close OMP_PLACES=cores
hostname
nvidia-smi --query-gpu=name,uuid --format=csv,noheader | tee gpu.txt
grep -q H200 gpu.txt
sha256sum -c inputs.sha256
for case_name in A1 B1 B2 A2; do
  TIMEFORMAT='%3R'; {{ time bash "$case_name/run.sh" 22; }} 2> "$case_name/elapsed_seconds.txt"
done
/home/kjhan/miniconda3/bin/python3 {audits[0]} {root}
'''
(root/'job.slurm').write_text(job);subprocess.run(['bash','-n',str(root/'job.slurm')],check=True)
print(root)

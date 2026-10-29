"""Stage/run a guarded continuation from a validated thermal64 t=.4125 checkpoint."""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_q1_checkpoint_geometry_20261029 import audit as geometry_audit
from audit_kepler_q1_thermal64_promotion_20261029 import payload, TAG
sys.path.insert(0,'/gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029')
from run_kepler_q1_orbital_cap_chain_20261029 import core_screen

root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64EarlyTo0825_H200_20261029')
parent=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Early_H200_20261029')
legacy=Path('/gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029')
ic=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/ic_v1/khout.000001.dat')
docs=Path('/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs')
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def save(p,v):
    with p.open('x') as f:json.dump(v,f,indent=2);f.write('\n')
def stage():
    pc=json.loads((parent/'contract.json').read_text())
    if pc['target_time']!=.4125 or pc['eta']!=.06 or sha(pc['binary'])!=pc['binary_sha256'] or sha(pc['source'])!=pc['source_sha256']:raise ValueError('Unexpected parent contract')
    if root.exists():raise FileExistsError(root)
    root.mkdir();shutil.copy2(__file__,root/'driver.py')
    files=[parent/'contract.json',parent/'run.sh',parent/'inputs.sha256',parent/'params.dat',Path(pc['binary']),Path(pc['source']),ic,legacy/'contract.json',legacy/'run_kepler_q1_orbital_cap_chain_20261029.py',docs/'audit_kepler_q1_thermal64_early_h200_20261029.py',docs/'audit_kepler_q1_bounded_gate_20261029.py',docs/'audit_kepler_q1_checkpoint_geometry_20261029.py',docs/'audit_kepler_q1_restart_dt_bundle_20261029.py',docs/'kepler_native_global_ids_20261029.py',docs/'audit_kepler_native_radial_20261026.py',docs/'audit_kepler_q1_thermal64_promotion_20261029.py',root/'driver.py']
    contract={'parent_job':439903,'parent':str(parent),'target':.825,'eta':.06,'binary_sha256':pc['binary_sha256'],
              'frozen':{str(p):sha(p) for p in files},
              'output_policy':'One ~52 MB final checkpoint and small audits; parent snapshots untouched.',
              'scope':'Gated continuation of accepted early native W0 thermal64 state; not t6.5 or long-time validation.'}
    save(root/'contract.json',contract)
    job='''#!/bin/bash
#SBATCH --job-name=KepThermEarly0825
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --dependency=afterok:439903
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --chdir='''+str(root)+'''
#SBATCH --output='''+str(root)+'''/slurm-%j.out
#SBATCH --error='''+str(root)+'''/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=8 OMP_STACKSIZE=256M OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
hostname
nvidia-smi --query-gpu=name,uuid --format=csv,noheader
export PYTHONPATH='''+str(docs)+'''
/home/kjhan/miniconda3/bin/python3 driver.py run
'''
    (root/'job.slurm').write_text(job)
    print(json.dumps(contract,indent=2))
def run():
    cfg=json.loads((root/'contract.json').read_text())
    if Path.cwd().resolve()!=root or not os.environ.get('SLURM_JOB_ID') or os.environ.get('SLURM_CPUS_PER_TASK')!='8':raise ValueError('Wrong allocation')
    for path,expected in cfg['frozen'].items():
        if sha(path)!=expected:raise ValueError('Changed frozen input '+path)
    if (parent/'exit_status.txt').read_text().strip()!='0':raise ValueError('Parent simulation failed')
    comparison=parent/'comparison.json'
    if not comparison.exists():subprocess.run([sys.executable,str(docs/'audit_kepler_q1_thermal64_early_h200_20261029.py')],check=True)
    audit=json.loads(comparison.read_text())
    pc=json.loads((parent/'contract.json').read_text());seedname=Path(pc['source']).name
    ends=[p for p in parent.glob('khout.*.dat') if p.name!=seedname]
    if len(ends)!=1 or not audit['passed'] or sha(ends[0])!=audit['candidate_endpoint_sha256'] or sha(parent/'eunha2.exe')!=cfg['binary_sha256']:raise ValueError('Parent checkpoint failed physical gates')
    endpoint=ends[0];p0,_=ordered_global(ic);p,meta=ordered_global(endpoint,p0['u4if']&ID_MASK);q=payload(endpoint,len(p))
    if meta['t']!=.4125 or not np.array_equal(p0['mass'],p['mass']) or np.any(p['w2']!=0) or np.any(q['tag']!=TAG) or np.any(~np.isfinite(q['u'])) or np.any(q['u']<=0):raise ValueError('Invalid thermal restart')
    step=int(re.fullmatch(r'khout\.(\d+)\.dat',endpoint.name)[1]);(root/endpoint.name).symlink_to(endpoint)
    (root/'eunha2.exe').symlink_to(Path(pc['binary']));shutil.copy2(parent/'params.dat',root/'params.dat')
    script=(parent/'run.sh').read_text().replace(str(parent),str(root))
    old='HYDRO_TSTOP=0.41249999999999998';new='HYDRO_TSTOP=0.825'
    if script.count(old)!=1:raise ValueError('Unexpected parent runner')
    script=script.replace(old,new);(root/'run.sh').write_text(script)
    inputs=[endpoint,Path(pc['binary']),root/'params.dat',root/'run.sh']
    (root/'inputs.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in inputs))
    save(root/'restart_receipt.json',{'source':str(endpoint),'source_sha256':sha(endpoint),'step':step,'parent_comparison_sha256':sha(comparison),'target':.825})
    subprocess.run(['bash','-n',str(root/'run.sh')],check=True)
    subprocess.run(['bash',str(root/'run.sh'),str(step)],cwd=root,check=True)
    final=[p for p in root.glob('khout.*.dat') if p.name!=endpoint.name]
    if len(final)!=1:raise ValueError('Ambiguous endpoint')
    gate_path=root/'checkpoint_gate.json'
    with (root/'checkpoint_gate.log').open('x') as f:
        gateproc=subprocess.run([sys.executable,str(docs/'audit_kepler_q1_bounded_gate_20261029.py'),str(ic),str(root),'--time','.825','--output',str(gate_path)],stdout=f,stderr=subprocess.STDOUT)
    gate=json.loads(gate_path.read_text())
    geom=geometry_audit(ic,final[0]);save(root/'independent_geometry.json',geom)
    core=core_screen(ic,final[0],geom,json.loads((legacy/'contract.json').read_text()));save(root/'core_screen.json',core)
    _,fm=ordered_global(final[0],p0['u4if']&ID_MASK);fq=payload(final[0],len(p0))
    checks={'native_exit_zero':(root/'exit_status.txt').read_text().strip()=='0','exact_time':fm['t']==.825,
            'original_28':gateproc.returncode==0 and len(gate['checks'])==28 and all(gate['checks'].values()),
            'geometry':geom['passed'],'core_screen':core['passed'],
            'thermal_payload':bool(np.all(fq['tag']==TAG) and np.all(np.isfinite(fq['u'])) and np.all(fq['u']>0) and np.all(np.isfinite(fq['v'])) and np.all(fq['v']>0))}
    result={'passed':all(checks.values()),'checks':checks,'endpoint':str(final[0]),'sha256':sha(final[0]),'time':fm['t'],
            'scope':'Gated early thermal64 continuation; not long-time stability.'}
    save(root/'result.json',result);print(json.dumps(result,indent=2))
    if not result['passed']:raise SystemExit(1)
if __name__=='__main__':
    if sys.argv[1]=='stage':stage()
    elif sys.argv[1]=='run':run()
    else:raise SystemExit('Expected stage or run')

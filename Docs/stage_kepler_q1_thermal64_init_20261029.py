"""Stage a guarded, zero-time thermal64 promotion diagnostic on bounded Q1."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate(path):
    cfg=json.loads(path.read_text());step=str(cfg['step'])
    expected={'GFS_THERMAL64':'1','GFS_THERMAL64_INIT_ONLY':'1','KEPLER_RESTART_STEP':step,
              'GFS_THERMAL64_INIT_TIME':'10','GFS_KEPLER_FIXED_POLYGON':'1',
              'GFS_KEPLER_WALL_EVENTS':'1','GFS_GPU_FINAL_VOLUME':'1','GFS_GPU_RAW_FACES':'1',
              'GFS_GPU_TESS_DUMP':'0','GFS_FAIL_NONPOSITIVE_U':'1','EUNHA_IC':'kepler',
              'SEDOV_PHASE1':'1','OMP_NUM_THREADS':'2','HYDRO_TSTOP':'10.0001'}
    for key,value in expected.items():
        if os.environ.get(key)!=value:raise ValueError('Unexpected promotion environment: '+key)
    forbidden=['GFS_RK4_AUDIT_STEP_MIN','GFS_RK4_AUDIT_STEP_MAX','GFS_RK4_AUDIT_IDS','GFS_RK4_AUDIT_STAGE',
               'GFS_GPU_FACE_DUMP_CALL','GFS_GPU_FACE_DUMP_IDS','GFS_GPU_FACE_DUMP_FIRST','GFS_GPU_FACE_DUMP_LAST',
               'GFS_RK4_DIRECT_U','GFS_RK4_TRANSACTION','GFS_DYN_W2','GFS_W2_LINEAR','GFS_ENTROPY_SWITCH',
               'GFS_RK4_SNAPSHOT_TEST_EXIT','GFS_VARPRESS_GEOM_AUDIT']
    if any(key in os.environ for key in forbidden):raise ValueError('Unrequested evolution or diagnostic mode')
    for key in ['initial','binary']:
        if sha(cfg[key])!=cfg[key+'_sha256']:raise ValueError('Frozen '+key+' changed')
    if sha(Path(__file__))!=cfg['stager_sha256']:raise ValueError('Promotion validator source changed')
    if sha(path.parent/'params.dat')!=cfg['params_sha256'] or sha(path.parent/'run.sh')!=cfg['runner_sha256']:
        raise ValueError('Promotion input or runner changed')


def stage(root):
    parent=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/chain_h200_v2')
    state=json.loads((parent/'chain_state.json').read_text());gate=json.loads((parent/'to10p0/checkpoint_gate.json').read_text())
    initial=Path(gate['endpoint']);binary=Path('/gpfs/kjhan/LagForce/Kepler_SlipWall_Integration_20261029/candidate_v2/eunha2.exe')
    if (state['status']!='finite_chain_complete' or gate['time']!=10. or len(gate['checks'])!=28
            or not all(gate['checks'].values()) or sha(initial)!=state['completed'][-1]['sha256']
            or sha(binary)!='f26816122df101ea71f72b4f3fb9233e965d06784ed94e88cfb46245fa251444'):
        raise ValueError('Requires exact accepted t10 native W0 checkpoint')
    if root.exists():raise FileExistsError('Refusing to overwrite promotion diagnostic')
    root.mkdir(parents=True);(root/initial.name).symlink_to(initial);(root/'eunha2.exe').symlink_to(binary)
    shutil.copy2(parent/'to10p0/params.dat',root/'params.dat');shutil.copy2(__file__,root/'promotion.py')
    step=int(initial.name.split('.')[1]);cfg={'initial':str(initial),'initial_sha256':sha(initial),'binary':str(binary),
        'binary_sha256':sha(binary),'step':step,'time':10.,'stager_sha256':sha(root/'promotion.py'),
        'output_policy':'One 52 MB tagged initialization snapshot; zero physical time advanced; no periodic or full-particle TSV.',
        'scope':'Initialization only. Thermal64 evolution and its wall-event/transaction compatibility remain explicitly unsupported.'}
    template=(parent/'to10p0/run.sh').read_text().replace('cd '+str(parent/'to10p0'),'cd '+str(root))
    template=template.replace('export OMP_NUM_THREADS=8','export OMP_NUM_THREADS=2').replace('export I_MPI_HYDRA_BOOTSTRAP=slurm','unset I_MPI_HYDRA_BOOTSTRAP')
    template=template.replace('HYDRO_TSTOP=10.0','HYDRO_TSTOP=10.0001')
    extra='unset GFS_RK4_DIRECT_U GFS_RK4_TRANSACTION GFS_RK4_SNAPSHOT_TEST_EXIT GFS_VARPRESS_GEOM_AUDIT GFS_GPU_FACE_DUMP_FIRST GFS_GPU_FACE_DUMP_LAST\n'
    extra+='export GFS_THERMAL64=1 GFS_THERMAL64_INIT_ONLY=1 KEPLER_RESTART_STEP='+str(step)+' GFS_THERMAL64_INIT_TIME=10\n'
    extra+='ulimit -f 131072\n/home/kjhan/miniconda3/bin/python3 '+str(root/'promotion.py')+' validate '+str(root/'promotion.json')+'\n'
    assert template.count('set +e\n')==1
    template=template.replace('set +e\n',extra+'set +e\n')
    (root/'run.sh').write_text(template)
    cfg.update(params_sha256=sha(root/'params.dat'),runner_sha256=sha(root/'run.sh'))
    (root/'promotion.json').write_text(json.dumps(cfg,indent=2)+'\n')
    job='''#!/bin/bash
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 CUDA_VISIBLE_DEVICES=GPU-42b3367f-1df7-2b2d-7961-7781c87b022d
bash '''+str(root/'run.sh')+' '+str(step)+'\n'
    (root/'job_ada.sh').write_text(job)
    print(json.dumps(cfg,indent=2))


if __name__=='__main__':
    if sys.argv[1]=='validate':validate(Path(sys.argv[2]))
    elif sys.argv[1]=='stage':stage(Path(sys.argv[2]).resolve())
    else:raise SystemExit('Expected stage or validate')

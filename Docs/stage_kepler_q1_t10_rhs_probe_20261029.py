"""Stage one bounded selected-ID native RHS probe, or validate its environment."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate(path):
    cfg=json.loads(path.read_text());step=str(cfg['restart_step'])
    expected={'GFS_WALL_RHS_PROBE_ONLY':'1','GFS_RK4_FAILURE_SNAPSHOT':'1',
              'GFS_RK4_AUDIT_STAGE':'k1_rhs','GFS_RK4_AUDIT_STEP_MIN':step,'GFS_RK4_AUDIT_STEP_MAX':step,
              'GFS_RK4_AUDIT_IDS':','.join(str(x['id']) for x in cfg['selected']),
              'GFS_KEPLER_FIXED_POLYGON':'1','GFS_KEPLER_WALL_EVENTS':'1',
              'GFS_GPU_TESS_DUMP':'0','GFS_GPU_RAW_FACES':'1','GFS_FAIL_NONPOSITIVE_U':'1',
              'EUNHA_KEPLER_DTETA':'0','EUNHA_IC':'kepler','SEDOV_PHASE1':'1',
              'HYDRO_TSTOP':'10.0001','OMP_NUM_THREADS':'2'}
    if len(cfg['selected'])!=8 or len({x['id'] for x in cfg['selected']})!=8:
        raise ValueError('Exactly eight distinct selected IDs required')
    for key,value in expected.items():
        if os.environ.get(key)!=value:raise ValueError('Unsafe diagnostic environment: '+key)
    forbidden=['GFS_GPU_FACE_DUMP_CALL','GFS_GPU_FACE_DUMP_IDS','GFS_GPU_FACE_DUMP_FIRST',
               'GFS_GPU_FACE_DUMP_LAST','GFS_RK4_SNAPSHOT_TEST_EXIT','GFS_VARPRESS_GEOM_AUDIT']
    if any(key in os.environ for key in forbidden):raise ValueError('Unrequested diagnostic enabled')
    if sha(cfg['initial'])!=cfg['initial_sha256'] or sha(cfg['endpoint'])!=cfg['endpoint_sha256']:
        raise ValueError('Accepted snapshot changed')
    if sha(cfg['binary'])!=cfg['binary_sha256']:raise ValueError('Native executable changed')
    for name,expected_hash in cfg['helper_hashes'].items():
        if sha(name)!=expected_hash:raise ValueError('Probe source changed: '+name)


def stage(root):
    base=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029');parent=base/'chain_h200_v2'
    initial=base/'ic_v1/khout.000001.dat';gate=json.loads((parent/'to10p0/checkpoint_gate.json').read_text())
    state=json.loads((parent/'chain_state.json').read_text());endpoint=Path(gate['endpoint'])
    if state['status']!='finite_chain_complete' or gate['time']!=10. or not all(gate['checks'].values()):
        raise ValueError('Requires the exact accepted t10 parent')
    if sha(endpoint)!=state['completed'][-1]['sha256']:raise ValueError('Accepted endpoint changed')
    a,_=ordered_global(initial);ids=a['u4if']&ID_MASK;b,meta=ordered_global(endpoint,ids)
    if meta['t']!=10. or not np.array_equal(a['mass'],b['mass']) or np.any(b['w2']!=0):
        raise ValueError('Unexpected checkpoint IDs, mass or W0')
    x=a['x']-3.;y=a['y']-3.;r=np.hypot(x,y);phase=10*(r*r+.01175**2)**(-.75)
    orbit=np.hypot(b['x']-(3+x*np.cos(phase)-y*np.sin(phase)),b['y']-(3+x*np.sin(phase)+y*np.cos(phase)))/np.maximum(r,.01175)
    density=np.abs(b['den'].astype(float)/a['den'].astype(float)-1)
    entropy=np.abs(np.log(b['ie'].astype(float)/a['ie'].astype(float))+(float(np.float32(1.666667))-1)*np.log(b['volume'].astype(float)/a['volume'].astype(float)))
    core=r<.1;adj=(r>=.1)&(r<.5);disk=(r>=.5)&(r<2.)
    specs=[('center',core,-r),('core_density_extreme',core,density),('core_orbit_extreme',core,orbit),
           ('core_entropy_extreme',core,entropy),('adjacent_density_extreme',adj,density),
           ('adjacent_orbit_extreme',adj,orbit),('disk_density_extreme',disk,density),('disk_entropy_extreme',disk,entropy)]
    chosen=[];selected=[]
    for label,mask,score in specs:
        candidates=np.flatnonzero(mask);rank=candidates[np.argsort(score[candidates])[::-1]]
        i=next(int(j) for j in rank if int(j) not in chosen);chosen.append(i)
        selected.append({'id':int(ids[i]),'purpose':label,'initial_radius':float(r[i]),
                         'density_change_abs':float(density[i]),'orbit_error_over_max_r0_eps':float(orbit[i])})
    binary=Path('/gpfs/kjhan/LagForce/Kepler_SlipWall_Integration_20261029/candidate_v2/eunha2.exe')
    if sha(binary)!='f26816122df101ea71f72b4f3fb9233e965d06784ed94e88cfb46245fa251444':raise ValueError('Native binary changed')
    if root.exists():raise FileExistsError('Refusing to overwrite probe')
    root.mkdir(parents=True);(root/endpoint.name).symlink_to(endpoint);(root/'eunha2.exe').symlink_to(binary)
    shutil.copy2(parent/'to10p0/params.dat',root/'params.dat');shutil.copy2(__file__,root/'probe.py')
    docs=Path(__file__).resolve().parent
    for name in ['kepler_native_global_ids_20261029.py','audit_kepler_native_radial_20261026.py']:
        shutil.copy2(docs/name,root/name)
    step=int(endpoint.name.split('.')[1]);cfg={'initial':str(initial),'initial_sha256':sha(initial),
        'endpoint':str(endpoint),'endpoint_sha256':sha(endpoint),'binary':str(binary),'binary_sha256':sha(binary),
        'restart_step':step,'selected':selected,'maximum_selected_rows':8,
        'helper_hashes':{str(p):sha(p) for p in [root/'probe.py',root/'kepler_native_global_ids_20261029.py',root/'audit_kepler_native_radial_20261026.py',Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/audit_snapshot.py'),Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/snapshot_offsets.json')]},
        'output_policy':'One bounded native K1 snapshot about 52 MB and eight selected-ID TSV rows. No accepted time evolution, no full-face dump. Per-process file-size ceiling at most 128 MiB.',
        'scope':'Single frozen-state Q1 native pressure acceleration/work sample. Purposefully includes extremes; not a population statistic, a force history, or proof that pressure is slow.'}
    (root/'probe.json').write_text(json.dumps(cfg,indent=2)+'\n')
    template=(base/'k1_preflight/run.sh').read_text()
    template=template.replace('cd '+str(base/'k1_preflight'),'cd '+str(root))
    template=template.replace('EUNHA_KEPLER_DTETA=0.06','EUNHA_KEPLER_DTETA=0')
    template=template.replace('HYDRO_TSTOP=1e-9','HYDRO_TSTOP=10.0001')
    template=template.replace('GFS_RK4_AUDIT_STEP_MIN=0 GFS_RK4_AUDIT_STEP_MAX=1000000',
                              'GFS_RK4_AUDIT_STEP_MIN='+str(step)+' GFS_RK4_AUDIT_STEP_MAX='+str(step))
    old='GFS_RK4_AUDIT_STAGE=k1_rhs GFS_RK4_AUDIT_IDS=74309,83540,94689,24770,69096,138444,78359,9036'
    template=template.replace(old,'GFS_RK4_AUDIT_STAGE=k1_rhs GFS_RK4_AUDIT_IDS='+','.join(str(x['id']) for x in selected))
    command='/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat 1 > log 2>&1'
    assert template.count(command)==1 and template.count('set +e\n')==1 and old not in template
    template=template.replace(command,command.replace('params.dat 1','params.dat '+str(step)))
    extra='unset GFS_GPU_FACE_DUMP_FIRST GFS_GPU_FACE_DUMP_LAST GFS_RK4_SNAPSHOT_TEST_EXIT GFS_VARPRESS_GEOM_AUDIT\n'
    extra+='ulimit -f 131072\n/home/kjhan/miniconda3/bin/python3 '+str(root/'probe.py')+' validate '+str(root/'probe.json')+'\n'
    template=template.replace('set +e\n',extra+'set +e\n')
    (root/'run.sh').write_text(template)
    job='''#!/bin/bash
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 CUDA_VISIBLE_DEVICES=GPU-42b3367f-1df7-2b2d-7961-7781c87b022d
bash '''+str(root/'run.sh')+'\n'
    (root/'job_ada.sh').write_text(job)
    print(json.dumps(cfg,indent=2))


if __name__=='__main__':
    if sys.argv[1]=='validate':validate(Path(sys.argv[2]))
    elif sys.argv[1]=='stage':stage(Path(sys.argv[2]).resolve())
    else:raise SystemExit('Expected stage or validate')

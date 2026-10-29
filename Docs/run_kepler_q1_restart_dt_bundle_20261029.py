"""One finite H200 qualification + timestep bundle after the accepted t=10 chain.

Never regenerates an IC, modifies the parent, or repairs inherited errors.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import time
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_const_config_cache_20261029 import audit as cache_audit


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)


def parent_checkpoint(cfg):
    parent=Path(cfg['parent_root'])
    state=json.loads((parent/'chain_state.json').read_text())
    if state.get('status')!='finite_chain_complete':
        raise ValueError('Parent chain is not complete; no t0 or earlier-time fallback')
    if state.get('job')!=cfg['parent_job']:raise ValueError('Different parent Slurm lineage')
    if digest(parent/'to10p0/params.dat')!=cfg['parent_params_sha256'] or digest(parent/'to10p0/eunha2.exe')!=cfg['binary_hashes'][cfg['baseline']]:
        raise ValueError('Parent physical parameters/executable differ from the qualified lineage')
    gate=json.loads((parent/'to10p0/checkpoint_gate.json').read_text())
    if gate.get('time')!=10. or len(gate.get('checks',{}))!=28 or not all(gate['checks'].values()):
        raise ValueError('Exact t10 and all original checkpoint gates are required')
    last=state['completed'][-1];endpoint=Path(gate['endpoint']).resolve()
    if last['time']!=10. or not last['all_checks_passed'] or Path(last['endpoint']).resolve()!=endpoint:
        raise ValueError('Parent accepted-checkpoint receipt disagrees')
    endpoint.relative_to((parent/'to10p0').resolve())
    if digest(endpoint)!=last['sha256']:raise ValueError('Accepted t10 snapshot changed')
    initial=Path(cfg['initial'])
    if digest(initial)!=cfg['initial_sha256']:raise ValueError('Scientific IC changed')
    p0,m0=ordered_global(initial);p,m=ordered_global(endpoint,p0['u4if']&ID_MASK)
    if (m['t']!=10. or m0['t']!=0. or len(p)!=90352
            or not np.array_equal(p['mass'],p0['mass']) or np.any(p['w2']!=0)
            or not all(np.isfinite(p[k]).all() and np.all(p[k]>0) for k in ['mass','volume','ie','pressure','den'])):
        raise ValueError('Accepted restart failed native ID/mass/state checks')
    match=re.fullmatch(r'khout\.(\d+)\.dat',endpoint.name)
    if not match or int(match[1])<=1:raise ValueError('Not a positive native continuation step')
    hs=[]
    for line in (parent/'to10p0/log').read_text(errors='replace').splitlines():
        if '[WALL_ACCEPT]'in line:hs.append(float(re.search(r'\bh=(\S+)',line)[1]))
    recent=np.array(hs[-65:-1])  # exclude the intentionally shortened endpoint
    if len(recent)<16 or not np.isfinite(recent).all() or np.any(recent<=0):
        raise ValueError('Need a complete recent physical timestep history')
    return {'endpoint':str(endpoint),'sha256':last['sha256'],'step':int(match[1]),
            'time':10.,'recent_h_min':float(recent.min()),'recent_h_median':float(np.median(recent)),
            'snapshot_bytes':endpoint.stat().st_size}


def validate_environment(expected):
    required={'GFS_KEPLER_WALL_EVENTS':'1','GFS_KEPLER_FIXED_POLYGON':'1',
              'GFS_GPU_RAW_FACES':'1','GFS_GPU_FINAL_VOLUME':'1','GFS_GPU_TESS_DUMP':'0',
              'GFS_FAIL_NONPOSITIVE_U':'1','SEDOV_PHASE1':'1','EUNHA_IC':'kepler',
              'OMP_NUM_THREADS':'8','GFS_EXACT_STOP':'1','GFS_RK4_GAP_PROFILE':'1',
              'GFS_NATIVE_CONST_CONFIG':str(expected['cache'])}
    for key,value in required.items():
        if os.environ.get(key)!=value:raise ValueError('Unexpected effective environment: '+key)
    forbidden=['GFS_RK4_AUDIT_STEP_MIN','GFS_RK4_AUDIT_STEP_MAX','GFS_RK4_AUDIT_STAGE','GFS_RK4_AUDIT_IDS',
               'GFS_GPU_FACE_DUMP_CALL','GFS_GPU_FACE_DUMP_IDS','GFS_GPU_FACE_DUMP_FIRST','GFS_GPU_FACE_DUMP_LAST',
               'GFS_THERMAL64','GFS_DYN_W2','GFS_W2_LINEAR','GFS_ENTROPY_SWITCH','GFS_DUAL_ENERGY',
               'GFS_HALF_LIMIT','GFS_WALL_RHS_PROBE_ONLY','GFS_DIAG_CHECKPOINT_STEP','GFS_DUMP_INITIAL']
    if any(key in os.environ for key in forbidden):raise ValueError('Unbounded diagnostics or incompatible physics enabled')
    if float(os.environ['HYDRO_TSTOP'])!=expected['target'] or float(os.environ['EUNHA_DUMP_DT'])!=1e9:
        raise ValueError('Unexpected output target/cadence')
    cap=expected['cap']
    if cap is None:
        if 'GFS_DIAG_DT_MAX'in os.environ:raise ValueError('Unrequested timestep ceiling')
    elif not (math.isfinite(cap) and cap>0 and float(os.environ.get('GFS_DIAG_DT_MAX','nan'))==cap):
        raise ValueError('Wrong timestep ceiling')


def native_case(root,cfg,seed,name,binary,cache,target,cap,state):
    stage=root/name
    if stage.exists():raise ValueError('Existing case output: explicit resume audit required')
    stage.mkdir(parents=True)
    for path,expected in cfg['binary_hashes'].items():
        if digest(path)!=expected:raise ValueError('Frozen executable changed')
    if digest(seed['endpoint'])!=seed['sha256']:raise ValueError('Restart changed')
    (stage/Path(seed['endpoint']).name).symlink_to(seed['endpoint'])
    (stage/'eunha2.exe').symlink_to(binary)
    (stage/'params.dat').write_bytes((root/'params.dat').read_bytes())
    expected={'target':target,'cache':cache,'cap':cap,'source':seed}
    save(stage/'case.json',expected)
    text=(root/'native_template.sh').read_text()
    text=text.replace('cd '+cfg['template_directory'],'cd '+shlex.quote(str(stage)))
    if 'HYDRO_TSTOP=10.0'not in text:raise ValueError('Unexpected frozen native template')
    text=text.replace('HYDRO_TSTOP=10.0','HYDRO_TSTOP='+format(target,'.17g'))
    extra='unset GFS_THERMAL64 GFS_DIAG_DT_MAX GFS_GPU_FACE_DUMP_FIRST GFS_GPU_FACE_DUMP_LAST\n'
    extra+='export GFS_RK4_GAP_PROFILE=1 GFS_NATIVE_CONST_CONFIG='+str(cache)+'\n'
    if cap is not None:extra+='export GFS_DIAG_DT_MAX='+format(cap,'.17g')+'\n'
    extra+=shlex.quote(sys.executable)+' '+shlex.quote(str(root/Path(__file__).name))+' validate-environment '+shlex.quote(str(stage/'case.json'))+'\n'
    if text.count('set +e\n')!=1:raise ValueError('Unexpected launch template')
    text=text.replace('set +e\n',extra+'set +e\n')
    command='/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat "$1" > log 2>&1'
    if text.count(command)!=1:raise ValueError('Unexpected native command')
    text=text.replace(command,'TIMEFORMAT="%R"\n{ time '+command+'; } 2> elapsed_seconds.txt')
    (stage/'run.sh').write_text(text)
    subprocess.run(['bash','-n',str(stage/'run.sh')],check=True)
    state.update(active_case=name,source=seed['endpoint']);save(root/'bundle_state.json',state)
    native=subprocess.Popen(['bash',str(stage/'run.sh'),str(seed['step'])],cwd=stage)
    state['launcher_pid']=native.pid;save(root/'bundle_state.json',state)
    rc=native.wait()
    if rc:raise RuntimeError('Native case failed; accepted parent and diagnostics preserved: '+name)
    outputs=[p for p in stage.glob('khout.*.dat') if p.name!=Path(seed['endpoint']).name]
    if len(outputs)!=1:raise ValueError('Missing unique native final checkpoint')
    _,meta=ordered_global(outputs[0])
    if meta['t']!=target:raise ValueError('Native exact target not reached')
    save(stage/'completed.json',{'endpoint':str(outputs[0]),'sha256':digest(outputs[0]),'time':meta['t']})
    state['completed'].append(name);save(root/'bundle_state.json',state)
    return outputs[0]


def run(root):
    root=root.resolve();cfg=json.loads((root/'bundle.json').read_text())
    if Path.cwd().resolve()!=root or not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_PER_TASK','0'))!=8:
        raise ValueError('Run only in the intended eight-CPU Slurm allocation')
    if os.environ['SLURM_JOB_ID']==cfg['parent_job']:raise ValueError('Do not run inside the still-owned parent allocation')
    for name,expected in cfg['frozen_hashes'].items():
        if digest(name)!=expected:raise ValueError('Frozen input/source changed: '+name)
    seed=parent_checkpoint(cfg)
    with (root/'execution_claim.json').open('x') as stream:json.dump({'job':os.environ['SLURM_JOB_ID'],'parent':seed},stream)
    state={'job':os.environ['SLURM_JOB_ID'],'parent_job':cfg['parent_job'],'status':'running','completed':[]}
    save(root/'bundle_state.json',state);save(root/'accepted_parent.json',seed);started=time.monotonic()
    try:
        qualification=root/'qualification';qualification.mkdir()
        target=10.+8*seed['recent_h_median']
        for name,binary,mode in [('A1',cfg['baseline'],0),('B1',cfg['candidate'],1),('B2',cfg['candidate'],1),('A2',cfg['baseline'],0)]:
            native_case(root,cfg,seed,'qualification/'+name,binary,mode,target,None,state)
        save(qualification/'stage.json',{'initial':seed['endpoint'],'initial_sha256':seed['sha256'],'target_time':target,
             'order':['A1','B1','B2','A2'],'baseline_binary_sha256':cfg['binary_hashes'][cfg['baseline']],
             'candidate_manifest':{'eunha2.exe':cfg['binary_hashes'][cfg['candidate']]},
             'predeclared_gates':{'ABBA_median_wall_per_step_speedup_min':1.05}})
        performance=cache_audit(qualification);save(qualification/'audit.json',performance)
        use_cache=performance['passed'];binary=cfg['candidate'] if use_cache else cfg['baseline']
        save(root/'selected_binary.json',{'cache_qualified_on_H200':use_cache,'binary':binary,
             'sha256':cfg['binary_hashes'][binary],'fallback_reason':None if use_cache else 'H200 parity/speed qualification failed; retain original v2'})
        wall=performance['ABBA_cached_median_wall'] if use_cache else performance['ABBA_baseline_median_wall']
        h=seed['recent_h_min'];duration=cfg['physical_interval']
        if not (math.isfinite(wall) and wall>0 and math.isfinite(duration) and duration>0):
            raise ValueError('Invalid timing or physical interval')
        estimate=2*wall*duration/h*7+120
        if estimate>cfg['allocation_seconds']-(time.monotonic()-started)-300:
            raise RuntimeError('Estimated timestep bundle exceeds remaining allocation; results retained, interval not shortened')
        controls=root/'controls';controls.mkdir()
        save(controls/'stage.json',{'parent':seed,'initial':cfg['initial'],'start_time':10.,'target_time':10.+duration,
             'requested_caps':{'coarse':None,'half':h/2,'quarter':h/4},'binary':binary,'cache':int(use_cache),
             'binary_sha256':cfg['binary_hashes'][binary],
             'scope':'Common accepted t10 restart. Measures incremental timestep sensitivity; inherited error is not reset or repaired.'})
        for name,cap in [('coarse',None),('half',h/2),('quarter',h/4)]:
            native_case(root,cfg,seed,'controls/'+name,binary,int(use_cache),10.+duration,cap,state)
            stage=controls/name
            with (stage/'checkpoint_gate.log').open('w') as out:
                result=subprocess.run([sys.executable,str(root/'audit_kepler_q1_bounded_gate_20261029.py'),cfg['initial'],str(stage),
                    '--time',str(10.+duration),'--output',str(stage/'checkpoint_gate.json')],stdout=out,stderr=subprocess.STDOUT)
            if result.returncode:raise RuntimeError('Original physical gate failed: '+name)
        with (controls/'comparison.log').open('w') as out:
            result=subprocess.run([sys.executable,str(root/'audit_kepler_q1_restart_dt_bundle_20261029.py'),str(controls),
                '--output',str(controls/'comparison.json')],stdout=out,stderr=subprocess.STDOUT)
        if result.returncode:raise RuntimeError('Timestep separation/finite-state evidence did not qualify; do not change thresholds')
        state['status']='finite_bundle_complete'
    except BaseException as exc:
        state.update(status='stopped',error=str(exc));save(root/'bundle_state.json',state);raise
    save(root/'bundle_state.json',state)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['run','inspect-parent','validate-environment']);parser.add_argument('path',type=Path)
    args=parser.parse_args()
    if args.mode=='validate-environment':validate_environment(json.loads(args.path.read_text()))
    elif args.mode=='inspect-parent':
        try:print(json.dumps(parent_checkpoint(json.loads((args.path/'bundle.json').read_text())),indent=2))
        except (ValueError,FileNotFoundError) as exc:print('NOT READY:',exc);raise SystemExit(3)
    else:run(args.path)

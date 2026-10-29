"""Native bounded Q1 orbital-cap comparison with checked restart checkpoints."""
import hashlib
import json
import os
from pathlib import Path
import shlex
import socket
import subprocess
import sys
import time
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_q1_checkpoint_geometry_20261029 import audit as geometry_audit
from audit_kepler_q1_restart_dt_bundle_20261029 import region_metrics
from audit_kepler_q1_gravity_cap_preflight_20261029 import screen


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)


def core_screen(initial, endpoint, geometry, contract):
    a,_=ordered_global(initial); ids=a['u4if']&ID_MASK
    b,meta=ordered_global(endpoint,ids)
    r=np.hypot(a['x']-3.,a['y']-3.);rows=[]
    for lo,hi in contract['regions_initial_radii']:
        mask=(r>=lo)&(r<hi);v=region_metrics(a,a,b,mask,meta['t'])
        geo=next(x for x in geometry['regions'] if x['rlo']==lo and x['rhi']==hi)
        values={'density_native_relative_to_circular_geometry':geo['measures']['native_density_relative_to_circular_geometry'],
                'circular_position_error_over_max_r0_eps':v['original_circular_position_error_over_max_r0_eps'],
                'log_entropy_relative_to_original_IC':v['original_t0_entropy_change']}
        checks={k:screen(values[k],limits) for k,limits in contract['core_accuracy_screen'].items()}
        rows.append({'rlo':lo,'rhi':hi,'cells':int(mask.sum()),'values':values,'checks':checks,
                     'passed':all(all(c.values()) for c in checks.values())})
    return {'passed':bool(geometry['passed'] and all(row['passed'] for row in rows)),
            'time':meta['t'],'regions':rows,
            'scope':'Preregistered numerical kinematic/entropy screen; not a complete GIZMO or boundary-convergence validation.'}


def run(root):
    root=root.resolve();cfg=json.loads((root/'chain.json').read_text())
    mode=cfg['execution_mode']
    if Path.cwd().resolve()!=root:raise ValueError('Wrong experiment working directory')
    if mode=='slurm':
        if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_PER_TASK','0'))!=8:
            raise ValueError('Requires the intended eight-CPU Slurm allocation')
        job=os.environ['SLURM_JOB_ID']
    elif mode=='ada':
        if (socket.gethostname().split('.')[0].lower()!='lageunha' or os.environ.get('SLURM_JOB_ID')
                or os.environ.get('OMP_NUM_THREADS')!='8'
                or os.environ.get('CUDA_VISIBLE_DEVICES')!=cfg['ada_gpu_uuid']):
            raise ValueError('Requires the explicitly authorized Ada host, GPU UUID and eight threads')
        job='ada-pid-'+str(os.getpid())
    else:raise ValueError('Unexpected native execution mode')
    for path,expected in cfg['frozen_hashes'].items():
        if sha(path)!=expected:raise ValueError('Frozen input changed: '+path)
    contract=json.loads((root/'contract.json').read_text())
    evidence=json.loads(Path(cfg['gravity_screen']).read_text())
    if (not evidence['passed'] or evidence['selected_eta']!=cfg['eta'] or
            evidence['restart_selection']['selected_source']!=cfg['initial'] or
            evidence['restart_selection']['selected_source_sha256']!=sha(cfg['initial'])):
        raise ValueError('The preregistered gravity/restart screens do not authorize this experiment')
    startup=json.loads(Path(cfg['startup_audit']).read_text())
    startup_core=json.loads(Path(cfg['startup_core_screen']).read_text())
    previous=Path(cfg['starting_checkpoint']);step=cfg['starting_step']
    a,_=ordered_global(cfg['initial']); ids=a['u4if']&ID_MASK
    b,meta=ordered_global(previous,ids)
    if (not startup['passed'] or not startup_core['passed'] or startup['orbital_eta']!=cfg['eta']
            or startup['endpoint']!=str(previous) or startup['endpoint_sha256']!=sha(previous)
            or startup['target_time']!=meta['t'] or startup_core['time']!=meta['t']
            or int(previous.name.split('.')[1])!=step or not np.array_equal(a['mass'],b['mass'])
            or np.any(b['w2']!=0)):
        raise ValueError('Accepted capped startup checkpoint or core validation failed')
    if cfg['cache_mode']:
        parity=json.loads(Path(cfg['cache_parity']).read_text())
        if not parity['passed'] or parity['candidate_binary_sha256']!=sha(cfg['binary']):
            raise ValueError('Eta-specific native cache parity did not pass')
    with (root/'execution_claim.json').open('x') as stream:
        json.dump({'job':job,'initial':cfg['initial'],
                   'starting_checkpoint':str(previous),'starting_sha256':sha(previous)},stream)
    state={'job':job,'status':'running','completed':[]}
    save(root/'chain_state.json',state)
    try:
        for stage in cfg['stages']:
            stage_dir=root/stage['directory'];stage_dir.mkdir()
            (stage_dir/previous.name).symlink_to(previous)
            (stage_dir/'eunha2.exe').symlink_to(cfg['binary'])
            (stage_dir/'params.dat').write_bytes((root/'params.dat').read_bytes())
            source_sha=sha(previous)
            save(stage_dir/'restart_receipt.json',{'source':str(previous),'sha256':source_sha,
                 'step':step,'target_time':stage['time'],'orbital_eta':cfg['eta']})
            expected={'target':stage['time'],'cache':cfg['cache_mode'],'cap':None}
            save(stage_dir/'case.json',expected)
            script=(root/'native_template.sh').read_text()
            script=script.replace('cd '+cfg['template_directory'],'cd '+shlex.quote(str(stage_dir)))
            script=script.replace('HYDRO_TSTOP=10.0','HYDRO_TSTOP='+format(stage['time'],'.17g'))
            script=script.replace('EUNHA_KEPLER_DTETA=0','EUNHA_KEPLER_DTETA='+format(cfg['eta'],'.17g'))
            if mode=='ada':script=script.replace('export I_MPI_HYDRA_BOOTSTRAP=slurm','unset I_MPI_HYDRA_BOOTSTRAP')
            if 'HYDRO_TSTOP=10.0'in script or 'EUNHA_KEPLER_DTETA=0\n'in script:
                raise ValueError('Native environment template did not render')
            extra='unset GFS_THERMAL64 GFS_DIAG_DT_MAX GFS_GPU_FACE_DUMP_FIRST GFS_GPU_FACE_DUMP_LAST\n'
            extra+='export GFS_RK4_GAP_PROFILE=1 GFS_NATIVE_CONST_CONFIG='+str(cfg['cache_mode'])+'\n'
            extra+='test "$EUNHA_KEPLER_DTETA" = '+shlex.quote(format(cfg['eta'],'.17g'))+'\n'
            extra+=shlex.quote(sys.executable)+' '+shlex.quote(cfg['environment_validator'])+' validate-environment '+shlex.quote(str(stage_dir/'case.json'))+'\n'
            if script.count('set +e\n')!=1:raise ValueError('Unexpected launch template')
            script=script.replace('set +e\n',extra+'set +e\n')
            (stage_dir/'run.sh').write_text(script)
            subprocess.run(['bash','-n',str(stage_dir/'run.sh')],check=True)
            state.update(active_stage=stage['directory'],restart=str(previous));save(root/'chain_state.json',state)
            started=time.monotonic()
            rc=subprocess.run(['bash',str(stage_dir/'run.sh'),str(step)],cwd=stage_dir).returncode
            if rc:raise RuntimeError('Native run failed; checkpoints preserved')
            outputs=[p for p in stage_dir.glob('khout.*.dat') if p.name!=previous.name]
            if len(outputs)!=1:raise ValueError('Native stage did not produce exactly one final checkpoint')
            endpoint=outputs[0]
            with (stage_dir/'checkpoint_gate.log').open('w') as stream:
                check=subprocess.run([sys.executable,str(root/'audit_kepler_q1_bounded_gate_20261029.py'),
                    cfg['initial'],str(stage_dir),'--time',str(stage['time']),
                    '--output',str(stage_dir/'checkpoint_gate.json')],stdout=stream,stderr=subprocess.STDOUT)
            gate=json.loads((stage_dir/'checkpoint_gate.json').read_text())
            if check.returncode or len(gate['checks'])!=28 or not all(gate['checks'].values()):
                raise RuntimeError('Original checkpoint gate failed')
            independent=geometry_audit(Path(cfg['initial']),endpoint)
            save(stage_dir/'independent_geometry.json',independent)
            core=core_screen(Path(cfg['initial']),endpoint,independent,contract)
            save(stage_dir/'core_screen.json',core)
            if not core['passed']:raise RuntimeError('Preregistered core screen failed')
            previous=endpoint;step=int(previous.name.split('.')[1])
            state['completed'].append({'time':stage['time'],'endpoint':str(previous),
                 'sha256':sha(previous),'restart_sha256':source_sha,'original_28_passed':True,
                 'core_screen_passed':True,'elapsed_seconds':time.monotonic()-started})
            save(root/'chain_state.json',state)
        state['status']='finite_chain_complete'
    except BaseException as exc:
        state.update(status='stopped',error=str(exc));save(root/'chain_state.json',state);raise
    save(root/'chain_state.json',state)


if __name__=='__main__':run(Path(sys.argv[1]))

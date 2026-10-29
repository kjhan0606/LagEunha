"""Check actual native orbital-cap enforcement and measured preflight cost."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(root):
    cfg=json.loads((root/'preflight.json').read_text())
    gate=json.loads((root/'checkpoint_gate.json').read_text())
    text=(root/'log').read_text(errors='replace');accept=[];times={}
    for line in text.splitlines():
        if '[WALL_ACCEPT]'in line:
            accept.append({'step':int(re.search(r'\bstep=(\d+)',line)[1]),
                           't':float(re.search(r't_start=(\S+)',line)[1]),
                           'h':float(re.search(r'\bh=(\S+)',line)[1])})
        if '[RK4 gap]'in line:
            times[int(re.search(r'\bstep=(\d+)',line)[1])]=float(re.search(r'wall=(\S+?)s',line)[1])
    h=cfg['universal_gravity_h_ceiling'];target=cfg['target_time']
    complete=bool(accept) and accept[0]['t']==0 and accept[-1]['t']+accept[-1]['h']==target
    complete &= all(abs(a['t']+a['h']-b['t'])<1e-15 for a,b in zip(accept[:-1],accept[1:]))
    regular=[a for a in accept if a['h']>.5*h]
    measured=[times[a['step']] for a in regular[2:]]
    if len(measured)<4:raise ValueError('Too few native timing samples')
    checks={'original_28_gates':len(gate['checks'])==28 and all(gate['checks'].values()),
            'exact_time':gate['time']==target,'complete_accepted_clock':bool(complete),
            'native_exit_zero':(root/'exit_status.txt').read_text().strip()=='0',
            'positive_h_under_orbital_ceiling':all(0<a['h']<=h*(1+1e-12) for a in accept),
            'twenty_full_orbital_cap_steps':len(regular)==20 and all(abs(a['h']/h-1)<1e-12 for a in regular),
            'frozen_initial':digest(Path(cfg['initial']))==cfg['initial_sha256'],
            'frozen_native':digest(root/'eunha2.exe')==cfg['native_binary_sha256'],
            'frozen_runner':digest(root/'run.sh')==cfg['run_sh_sha256'],
            'no_full_particle_TSV':not(root/'rk4_stage_r0.tsv').exists()}
    wall=float(np.median(measured)); elapsed=float((root/'elapsed_seconds.txt').read_text().strip())
    return {'passed':all(checks.values()),'checks':{k:bool(v) for k,v in checks.items()},
            'accepted_steps':len(accept),'regular_steps':len(regular),'roundoff_endpoint_tails':[a for a in accept if a['h']<=.5*h],
            'target_time':target,'orbital_eta':cfg['orbital_eta'],'maximum_accepted_h':max(a['h'] for a in accept),
            'native_RK4_median_wall_seconds':wall,'timing_samples_after_two_warmup_steps':len(measured),
            'native_exec_wall_seconds':elapsed,'endpoint':gate['endpoint'],'endpoint_sha256':digest(Path(gate['endpoint'])),
            'frozen_cost_extrapolation':{'minimum_cap_steps_to_t10':int(np.ceil(10/h)),
                'same_median_step_cost_hours_to_t10':float(np.ceil(10/h)*wall/3600),
                'assumptions':'Ada original v2, frozen early geometry/per-step cost, cap dominant. Not a measured long run, H200 estimate, upper bound, or claimed cached speedup.'},
            'source_hashes':{str(p):digest(p) for p in [root/'preflight.json',root/'checkpoint_gate.json',root/'log',Path(__file__)]},
            'scope':'Actual native CUDA+MPI+OpenMP orbital-cap setup/timing check only, about 20 steps. One final IEEE-roundoff time tail is reported explicitly. Does not establish cap accuracy through .825/3/6.5/10, final physical stability, or GIZMO runtime competitiveness.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=audit(args.root)
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))
    if not result['passed']:raise SystemExit(1)

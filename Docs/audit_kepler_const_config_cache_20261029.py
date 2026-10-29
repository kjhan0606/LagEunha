"""Full-record/ledger parity and ABBA native constant-config cache timing."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from audit_kepler_restart_bytes_20261029 import compare
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK


def audit(root):
    stage=json.loads((root/'stage.json').read_text())
    initial,m0=ordered_global(Path(stage['initial']));ids=initial['u4if']&ID_MASK
    endpoints={};timings={};ledgers={};checks={}
    checks['frozen_initial_hash']=hashlib.sha256(Path(stage['initial']).read_bytes()).hexdigest()==stage['initial_sha256']
    for name in stage['order']:
        d=root/name;text=(d/'log').read_text(errors='replace')
        files=[p for p in d.glob('khout.*.dat') if p.name!=Path(stage['initial']).name]
        if len(files)!=1:raise ValueError('Expected one complete endpoint: '+name)
        p,m=ordered_global(files[0],ids);endpoints[name]=files[0]
        expected_binary=(stage.get('baseline_binary_sha256','f26816122df101ea71f72b4f3fb9233e965d06784ed94e88cfb46245fa251444')
                         if name in ['A1','A2'] else stage['candidate_manifest']['eunha2.exe'])
        checks[name+'_frozen_binary']=hashlib.sha256((d/'eunha2.exe').read_bytes()).hexdigest()==expected_binary
        ledger=[s for s in text.splitlines() if '[RK4E]'in s or '[WALL_ACCEPT]'in s]
        ledgers[name]=ledger
        rows=[]
        for line in text.splitlines():
            if '[RK4 gap]'not in line:continue
            rows.append({key:float(re.search(r'\b'+key+r'=(\S+?)s(?:\s|$)',line)[1]) for key in ['wall','tracked','prelude','after_update','after_force','after_post','stage_diag']})
        if len(rows)<4:raise ValueError('Too few complete native timing steps')
        measured=rows[2:]
        elapsed=float((d/'elapsed_seconds.txt').read_text().strip())
        timings[name]={'accepted_steps':len(rows),'warmup_steps_excluded':2,
                       'median_seconds':{key:float(np.median([r[key] for r in measured])) for key in measured[0]},
                       'native_exec_wall_seconds':elapsed,'native_exec_wall_per_physical_time':elapsed/(m['t']-m0['t'])}
        zero=all(float(re.search(r'\b'+key+r'=(\S+)',line)[1])==0 for line in ledger if '[RK4E]'in line for key in ['floor_cum','de_cum','sfl_cum'])
        checks[name+'_native_state_and_policy']=((d/'exit_status.txt').read_text().strip()=='0' and m['t']==stage['target_time']
            and np.array_equal(initial['mass'],p['mass']) and np.all(p['w2']==0) and zero
            and all(np.isfinite(p[k]).all() and np.all(p[k]>0) for k in ['mass','volume','ie','den','pressure'])
            and all(np.isfinite(p[k]).all() for k in ['x','y','vx','vy'])
            and not (d/'rk4_stage_r0.tsv').exists() and '[WALL_EVENT_STOP]'not in text)
    parity={name:compare(endpoints['A1'],endpoints[name]) for name in stage['order'][1:]}
    checks['all_full_restart_parities']=all(p['passed'] for p in parity.values())
    checks['all_complete_numeric_ledgers_exact']=all(ledgers[name]==ledgers['A1'] for name in stage['order'][1:])
    baseline=float(np.median([timings[n]['median_seconds']['wall'] for n in ['A1','A2']]))
    cached=float(np.median([timings[n]['median_seconds']['wall'] for n in ['B1','B2']]))
    speedup=baseline/cached
    checks['ABBA_median_wall_speedup_gate']=speedup>=stage['predeclared_gates']['ABBA_median_wall_per_step_speedup_min']
    return {'passed':all(checks.values()),'checks':{k:bool(v) for k,v in checks.items()},'full_record_parity':parity,
            'timings':timings,'ABBA_baseline_median_wall':baseline,'ABBA_cached_median_wall':cached,'ABBA_wall_speedup':speedup,
            'initial_sha256':hashlib.sha256(Path(stage['initial']).read_bytes()).hexdigest(),
            'scope':'Actual native restart ABBA comparison with identical physical time and steps. Excludes two startup steps from per-step timing; full native execution timing is also reported. Does not prove scaling, thermal64-enabled performance, long-time stability or performance at different geometry.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=audit(args.root)
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))
    if not result['passed']:raise SystemExit(1)

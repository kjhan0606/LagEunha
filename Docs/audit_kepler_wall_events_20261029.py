"""Audit a native late-state wall-event replay, without claiming disk stability."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from audit_kepler_native_radial_20261026 import ordered,energy


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(root, reader=ordered):
    stage=json.loads((root/'stage.json').read_text())
    initial=root/'khout.011580.dat';a,m0=reader(initial)
    endpoint=None
    for p in root.glob('khout.*.dat'):
        if p.name==initial.name:continue
        _,meta=reader(p)
        if meta['t']==stage['target_time']:
            if endpoint is not None:raise ValueError('Repeated exact-time endpoint')
            endpoint=p
    if endpoint is None:raise ValueError('Missing complete exact-time endpoint')
    b,meta=reader(endpoint);text=(root/'log').read_text(errors='replace')
    impacts=[]
    pattern=r'\[WALL_IMPACT\] step=(\d+) t=(\S+) id=(\d+) side=(\d+) gap=(\S+) tolerance=(\S+) vn_before=(\S+) kinetic_roundoff=(\S+)'
    for q in re.finditer(pattern,text):
        impacts.append({'step':int(q[1]),'time':float(q[2]),'id':int(q[3]),'side':int(q[4]),
                        'gap':float(q[5]),'tolerance':float(q[6]),'vn_before':float(q[7]),'kinetic_roundoff':float(q[8])})
    accepts=[q for q in text.splitlines() if '[WALL_ACCEPT]' in q]
    moment=[];angular=[];durations=[];starts=[]
    for q in accepts:
        m=re.search(r'momentum_residual=\(([^,]+),([^\)]+)\)',q)
        if not m:raise ValueError('Missing wall momentum ledger')
        moment.append(float(np.hypot(float(m[1]),float(m[2]))))
        angular.append(abs(float(re.search(r'angular_residual=(\S+)',q)[1])))
        durations.append(float(re.search(r'\bh=(\S+)',q)[1]))
        starts.append(float(re.search(r'\bt_start=(\S+)',q)[1]))
    ledger=[q for q in text.splitlines() if '[RK4E]' in q]
    injection_free=all(all(float(re.search(r'\b'+k+r'=(\S+)',q)[1])==0.
                          for k in ['floor_cum','de_cum','sfl_cum']) for q in ledger)
    x=b['x'].astype(float)-3.;y=b['y'].astype(float)-3.
    theta=np.arctan2(y,x);nearest=np.rint(theta/(2*np.pi/1024)-.5).astype(int)
    apothem=2.75*np.cos(np.pi/1024)
    gap=np.full(len(b),-np.inf)
    for shift in [-1,0,1]:
        angle=(nearest+shift+.5)*2*np.pi/1024
        gap=np.maximum(gap,x*np.cos(angle)+y*np.sin(angle)-apothem)
    tolerance=64*np.finfo(float).eps*np.maximum(2.75,abs(x)+abs(y))
    e0=energy(a,'baseline');et=energy(b,'baseline');de=(et-e0)/abs(e0)
    gates=stage['predeclared_gates']
    checks={'native_exit_zero':(root/'exit_status.txt').read_text().strip()=='0',
            'native_binary_and_initial_hashes':digest(root/'eunha2.exe')==stage['binary_sha256'] and digest(initial)==stage['input_sha256'],
            'complete_owned_IDs':len(a)==len(b)==55351,
            'same_fixed_masses':bool(np.array_equal(a['mass'],b['mass'])),
            'W0_exact':bool(np.all(a['w2']==0) and np.all(b['w2']==0)),
            'finite_positive_state':all(np.isfinite(b[k]).all() and np.all(b[k]>0) for k in ['mass','volume','ie','pressure','den']),
            'all_endpoints_inside_polygon_to_roundoff':bool(np.all(gap<=tolerance)),
            'no_energy_injection':bool(ledger) and injection_free,
            'expected_collision_reached':any((q['id'],q['side'])==tuple(gates['collision_ID_side']) for q in impacts),
            'all_impact_gaps_within_tolerance':all(abs(q['gap'])<=q['tolerance'] and q['vn_before']>0 for q in impacts),
            'one_permanent_ledger_per_accepted_prefix':len(ledger)==len(accepts)>0,
            'positive_physical_time_prefixes':all(h>0 for h in durations) and all(starts[i+1]>starts[i] for i in range(len(starts)-1)),
            'prefix_times_close':bool(starts) and abs((starts[-1]+durations[-1])-meta['t'])<1.e-12 and starts[0]==m0['t'],
            'wall_momentum_and_torque_balance':bool(moment) and max(moment+angular)<gates['max_step_momentum_and_angular_residual'],
            'restart_total_energy':abs(de)<gates['restart_relative_total_energy_change_max'],
            'no_failed_rollback':not re.search(r'payload_counts_dtold_ledgers_exact=0|device_cache_invalidated=0',text),
            'no_stop_or_original_crossing_abort':'[WALL_EVENT_STOP]' not in text and '[FIXED_POLYGON_CROSSING]' not in text}
    return {'initial':str(initial),'initial_sha256':digest(initial),'initial_time':m0['t'],
            'endpoint':str(endpoint),'endpoint_sha256':digest(endpoint),'time':meta['t'],
            'checks':{k:bool(v) for k,v in checks.items()},'passed':all(checks.values()),
            'impacts':impacts,'accepted_prefixes':len(accepts),'rollbacks':text.count('[WALL_ROLLBACK]'),
            'minimum_accepted_h':min(durations) if durations else None,
            'maximum_momentum_residual':max(moment,default=None),'maximum_angular_residual':max(angular,default=None),
            'restart_relative_total_energy_change':de,'minimum_internal_energy':float(b['ie'].min()),
            'maximum_endpoint_plane_penetration':float(gap.max()),
            'scope':'Late bounded Q2 event replay past original failure. Validates local collision/rollback/ledger behavior only; not Q1 boundary IC, full event-order convergence, long-time stability or GIZMO reproduction.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=analyze(args.root)
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(r,indent=2))
    if not r['passed']:raise SystemExit(1)

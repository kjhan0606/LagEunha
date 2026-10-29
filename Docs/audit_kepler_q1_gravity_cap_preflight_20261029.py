"""Restart eligibility and gravity-only cap screen on saved bounded Q1 data.

No hydrodynamic evolution, IC regeneration or native-state modification.
"""
import argparse
import hashlib
import json
from pathlib import Path
import socket
import time
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_q1_restart_dt_bundle_20261029 import region_metrics
from audit_kepler_q1_core_rk4_control_20261029 import advance
from audit_kepler_native_radial_20261026 import mean

EPS = .01175


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def screen(values, limits):
    return {key: bool(values[key.removesuffix('_max')] < limit)
            for key, limit in limits.items()}


def run(root):
    if socket.gethostname().split('.')[0].lower() in ('grammar', 'syntax'):
        raise RuntimeError('Run on the compute host')
    started = time.monotonic(); contract = json.loads((root/'contract.json').read_text())
    initial = Path(contract['initial']); candidate = Path(contract['earliest_chain_restart_candidate'])
    a, ma = ordered_global(initial); ids = a['u4if'] & ID_MASK
    b, mb = ordered_global(candidate, ids)
    gate = json.loads((candidate.parent/'checkpoint_gate.json').read_text())
    geometry = json.loads((root/'geometry_t0825.json').read_text())
    if (ma['t'] != 0 or mb['t'] != .825 or len(gate['checks']) != 28
            or not all(gate['checks'].values()) or not geometry['passed']
            or geometry['time'] != mb['t'] or geometry['source_hashes'][str(candidate)] != digest(candidate)
            or geometry['source_hashes'][str(initial)] != digest(initial)
            or not np.array_equal(a['mass'],b['mass'])):
        raise ValueError('Native checkpoint/source/independent geometry lineage failed')
    r0 = np.hypot(a['x']-3.,a['y']-3.); regions = []; eligible = True
    for lo, hi in contract['regions_initial_radii']:
        mask = (r0>=lo)&(r0<hi)
        actual = region_metrics(a,a,b,mask,mb['t'])
        geom = next(v for v in geometry['regions'] if v['rlo']==lo and v['rhi']==hi)
        values = {'density_native_relative_to_circular_geometry':geom['measures']['native_density_relative_to_circular_geometry'],
                  'circular_position_error_over_max_r0_eps':actual['original_circular_position_error_over_max_r0_eps'],
                  'log_entropy_relative_to_original_IC':actual['original_t0_entropy_change']}
        checks = {name:screen(values[name],limits) for name,limits in contract['core_accuracy_screen'].items()}
        passed = all(all(v.values()) for v in checks.values()); eligible &= passed
        regions.append({'rlo':lo,'rhi':hi,'cells':int(mask.sum()),'passed':passed,'checks':checks,'values':values})
    selected = candidate if eligible else initial
    selection = {'passed':True,'candidate_eligible':bool(eligible),'candidate':str(candidate),'regions':regions,
                 'selected_source':str(selected),'selected_source_sha256':digest(selected),
                 'selected_time':mb['t'] if eligible else 0.,
                 'reason':'Earliest chain checkpoint satisfies preregistered core screens' if eligible else
                          'Earliest chain checkpoint fails added core screens; reuse exact saved original IC, never regenerate or select a later state to hide failure'}
    with (root/'restart_eligibility.json').open('x') as stream:
        json.dump(selection,stream,indent=2);stream.write('\n')
    # Only independent central-gravity ODEs, not the native gas system.
    selected_core = r0<.5; radius = r0[selected_core]; mass = a['mass'][selected_core].astype(float)
    x = a['x'][selected_core]-3.;y = a['y'][selected_core]-3.
    q0 = np.column_stack((x,y,a['vx'][selected_core],a['vy'][selected_core]))
    center = radius==0
    if not np.any(center) or np.any(q0[center]!=0):
        raise ValueError('Expected stationary central generator in the gravity-only reference')
    omega = (radius*radius+EPS*EPS)**(-.75); cases=[]
    budget = contract['gravity_only_screen']['position_error_budget_fraction_of_native_screen']
    limits = {k:v*budget for k,v in contract['core_accuracy_screen']['circular_position_error_over_max_r0_eps'].items()}
    for eta in contract['candidate_eta_descending']:
        q = q0.copy(); t = 0.; hmax = eta*EPS**1.5; rows=[]; steps=0
        for target in contract['gravity_only_screen']['targets']:
            count = int(np.floor((target-t)/hmax))
            for _ in range(count):q = advance(q,hmax)
            tail = target-(t+count*hmax)
            if tail>0:q = advance(q,tail)
            steps += count+int(tail>0);t=target
            phase = t*omega
            exact = np.column_stack((x*np.cos(phase)-y*np.sin(phase),x*np.sin(phase)+y*np.cos(phase)))
            error = np.linalg.norm(q[:,:2]-exact,axis=1)/np.maximum(radius,EPS)
            details=[]
            for lo, hi in contract['regions_initial_radii']:
                mask=(radius>=lo)&(radius<hi);v=error[mask];m=mass[mask]
                values={'absolute_mass_mean':mean(v,m),'absolute_cell_p99':float(np.percentile(v,99)),
                        'absolute_max':float(np.max(v))}
                checks=screen(values,limits)
                details.append({'rlo':lo,'rhi':hi,'values':values,'checks':checks,'passed':all(checks.values())})
            rows.append({'time':target,'steps':steps,'finite':bool(np.isfinite(q).all()),'regions':details,
                         'passed':bool(np.isfinite(q).all()) and all(v['passed'] for v in details)})
            print('gravity-only eta=',eta,'time=',target,'steps=',steps,'passed=',rows[-1]['passed'],flush=True)
        cases.append({'eta':eta,'global_h_ceiling':hmax,'targets':rows,'passed':all(v['passed'] for v in rows)})
    qualified = [v for v in cases if v['passed']]
    eta = max(v['eta'] for v in qualified) if qualified else None
    result = {'passed':eta is not None,'restart_selection':selection,'selected_eta':eta,'cases':cases,
              'contract_sha256':digest(root/'contract.json'),
              'source_hashes':{str(p):digest(p) for p in [initial,candidate,root/'geometry_t0825.json',Path(__file__),
                  Path(__file__).with_name('audit_kepler_q1_core_rk4_control_20261029.py'),
                  Path(__file__).with_name('audit_kepler_q1_restart_dt_bundle_20261029.py')]},
              'analysis_seconds':time.monotonic()-started,
              'scope':'Gravity-only RK4 ODE screening from the preserved original IC on all original r<.5 sites. A stationary central generator makes eta*eps^1.5 the global orbital ceiling. This supports a candidate native cap through t10 only. It neither evolves gas, repairs inherited error, proves long-time accuracy, nor qualifies t120/600. Current native chain and dependent controls are unchanged.'}
    with (root/'gravity_cap_screen.json').open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    args=parser.parse_args();result=run(args.root)
    print('Candidate native eta:',result['selected_eta'],'source:',result['restart_selection']['selected_source'],flush=True)
    if not result['passed']:raise SystemExit(1)

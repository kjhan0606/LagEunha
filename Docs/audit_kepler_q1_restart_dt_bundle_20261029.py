"""Three native timestep controls from one inherited accepted Q1 checkpoint."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_native_radial_20261026 import mean, energy, GAMMA_MINUS_ONE


def region_metrics(initial, seed, final, mask, target_time):
    mass=seed['mass'][mask].astype(float)
    def describe(values):
        a=np.asarray(values)[mask]
        return {'signed_mass_mean':mean(a,mass),'absolute_mass_mean':mean(abs(a),mass),
                'absolute_cell_p99':float(np.percentile(abs(a),99)),
                'absolute_max':float(np.max(abs(a)))}
    def angular(state):
        return float(np.sum(mass*((state['x'][mask]-3.)*state['vy'][mask]
                            -(state['y'][mask]-3.)*state['vx'][mask]),dtype=np.longdouble))
    entropy=np.log(final['ie'].astype(float)/seed['ie'].astype(float))+GAMMA_MINUS_ONE*np.log(final['volume'].astype(float)/seed['volume'].astype(float))
    entropy0=np.log(final['ie'].astype(float)/initial['ie'].astype(float))+GAMMA_MINUS_ONE*np.log(final['volume'].astype(float)/initial['volume'].astype(float))
    x0=initial['x']-3.;y0=initial['y']-3.;r0=np.hypot(x0,y0)
    phase=target_time*(r0*r0+.01175**2)**(-.75)
    position=np.hypot(final['x']-(3.+x0*np.cos(phase)-y0*np.sin(phase)),
                      final['y']-(3.+x0*np.sin(phase)+y0*np.cos(phase)))
    e0=energy(seed[mask],'baseline');l0=angular(seed)
    if not np.isfinite(e0) or e0==0 or not np.isfinite(l0) or l0==0:
        raise ValueError('Undefined regional invariant normalization')
    return {'cells':int(mask.sum()),'incremental_entropy_signed':mean(entropy[mask],mass),
            'incremental_entropy_abs':mean(abs(entropy[mask]),mass),
            'incremental_density_abs':mean(abs(final['den'][mask].astype(float)/seed['den'][mask].astype(float)-1),mass),
            'incremental_entropy_distribution':describe(entropy),
            'restart_relative_energy_change':(energy(final[mask],'baseline')-e0)/abs(e0),
            'restart_relative_angular_momentum_change':angular(final)/l0-1,
            'original_t0_density_change':describe(final['den'].astype(float)/initial['den'].astype(float)-1),
            'original_t0_entropy_change':describe(entropy0),
            'original_circular_position_error':describe(position),
            'original_circular_position_error_over_max_r0_eps':describe(position/np.maximum(r0,.01175))}


def audit(root):
    cfg=json.loads((root/'stage.json').read_text())
    seed=Path(cfg['parent']['endpoint']);p0,m0=ordered_global(seed);ids=p0['u4if']&ID_MASK
    if hashlib.sha256(seed.read_bytes()).hexdigest()!=cfg['parent']['sha256']:
        raise ValueError('Accepted parent changed')
    if m0['t']!=cfg['start_time'] or m0['t']!=10.:
        raise ValueError('Requires the exact accepted t10 restart')
    initial,_=ordered_global(Path(cfg['initial']),ids)
    r0=np.hypot(initial['x']-3.,initial['y']-3.)
    masks={'core_r0_lt_0p1':r0<.1,'core_r0_0p1_to0p5':(r0>=.1)&(r0<.5),
           'disk_r0_0p5_to2':(r0>=.5)&(r0<2.),'exterior_r0_ge2':r0>=2.}
    states={};cases={};checks={}
    for name in ['coarse','half','quarter']:
        d=root/name;gate=json.loads((d/'checkpoint_gate.json').read_text())
        checks[name+'_same_restart_bytes']=hashlib.sha256((d/seed.name).read_bytes()).hexdigest()==cfg['parent']['sha256']
        checks[name+'_same_frozen_executable']=hashlib.sha256((d/'eunha2.exe').read_bytes()).hexdigest()==cfg['binary_sha256']
        p,m=ordered_global(Path(gate['endpoint']),ids);states[name]=p
        text=(d/'log').read_text(errors='replace');rows=[]
        for line in text.splitlines():
            if '[WALL_ACCEPT]'in line:
                rows.append((float(re.search(r't_start=(\S+)',line)[1]),float(re.search(r'\bh=(\S+)',line)[1])))
        hs=np.array([v[1] for v in rows]);cap=cfg['requested_caps'][name]
        covered=(bool(rows) and rows[0][0]==m0['t'] and abs(rows[-1][0]+rows[-1][1]-m['t'])<1e-12
                 and all(abs(t+h-rows[i+1][0])<1e-12 for i,(t,h) in enumerate(rows[:-1])))
        checks[name+'_all_original_physics_gates']=len(gate['checks'])==28 and all(gate['checks'].values())
        checks[name+'_same_parent_and_exact_target']=covered and m['t']==cfg['target_time'] and np.array_equal(p0['mass'],p['mass'])
        checks[name+'_positive_h_and_correct_cap']=bool(len(hs)>0 and np.all(np.isfinite(hs)) and np.all(hs>0) and (cap is None or np.all(hs<=cap*(1+1e-12))))
        regions={}
        for label,mask in masks.items():
            regions[label]=region_metrics(initial,p0,p,mask,cfg['target_time'])
        cases[name]={'endpoint':gate['endpoint'],'sha256':hashlib.sha256(Path(gate['endpoint']).read_bytes()).hexdigest(),
            'accepted_steps':len(hs),'h_min':float(hs.min()),'h_median':float(np.median(hs)),
            'h_max':float(hs.max()),'requested_cap':cap,'restart_relative_energy_change':(energy(p,'baseline')-energy(p0,'baseline'))/abs(energy(p0,'baseline')),
            'initial_t0_absolute_gate_metrics':{key:gate[key] for key in ['disk_density_abs_mass_weighted','disk_entropy_log_abs_mass_weighted','disk_position_error_abs_mass_weighted','total_energy_relative']},
            'restart_increment_regions':regions,'native_exec_wall_seconds':float((d/'elapsed_seconds.txt').read_text().strip())}
    ratios={'half_over_coarse':cases['half']['h_median']/cases['coarse']['h_median'],
            'quarter_over_half':cases['quarter']['h_median']/cases['half']['h_median']}
    checks['actual_timestep_levels_separated']=ratios['half_over_coarse']<.75 and ratios['quarter_over_half']<.75
    comparisons={}
    for name in ['coarse','half']:
        fields={}
        for label,mask in masks.items():
            mass=p0['mass'][mask].astype(float);row={}
            for key in ['x','y','vx','vy','volume','ie','den','pressure']:
                a=states[name][key][mask].astype(float);b=states['quarter'][key][mask].astype(float)
                row[key]={'mass_RMS_absolute':float(np.sqrt(mean((a-b)**2,mass))),
                          'max_absolute':float(np.max(abs(a-b)))}
            fields[label]=row
        comparisons[name+'_versus_quarter']=fields
    return {'passed':bool(all(checks.values())),'checks':{k:bool(v) for k,v in checks.items()},'cases':cases,
        'observed_median_h_ratios':ratios,'endpoint_differences':comparisons,'start_time':m0['t'],'time':cfg['target_time'],
        'scope':'Three real native trajectories over the same short physical interval, from the identical inherited accepted t10 state. Reports incremental sensitivity and absolute old-IC metrics. Does not reset inherited core error, establish formal convergence order, prove improved thermal accuracy, or qualify a t120/600 expansion by itself.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=audit(args.root)
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))
    if not result['passed']:raise SystemExit(1)

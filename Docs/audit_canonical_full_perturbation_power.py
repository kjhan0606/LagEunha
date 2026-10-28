"""Close pressure+TCF+HG+shear tendency against completed RK4 split steps."""
import importlib.util
import json
from pathlib import Path

import numpy as np
from audit_canonical_pair_history import gravity_coefficient

base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_pair_hg_modes_20261018'
spec=importlib.util.spec_from_file_location('reader',base/'run_canonical_av_branches_20261018.py')
reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
activation=json.loads((base/'Codex_pair_history_20261018/activation_v2.json').read_text())
power=activation['cases']['t6']['actual']['actual_velocity_perturbation_power']
hg=json.loads((root/'analysis.json').read_text())['states']['actual']['regions']
start,t0=reader.read(root/'khout.000656.dat')
initial,_=reader.read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
r0=np.hypot(initial['x']-3,initial['y']-3)
masks={'initial_inner':(r0<.5)&(initial['den']<.5),
       'initial_inner_edge':(r0<.6)&(initial['den']>.5),
       'initial_disk_core':(r0>=.6)&(r0<1.9)&(initial['den']>.5)}
def cell_energy(a):
    x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']))
    omega=np.sqrt(gravity_coefficient(np.linalg.norm(x,axis=1)))
    delta=v-omega[:,None]*np.column_stack((-x[:,1],x[:,0]))
    return .5*a['mass']*np.sum(delta**2,axis=1)
energy0=cell_energy(start)
out={'scope':'Instantaneous complete perturbation kinetic-power budget at t6 plus finite-dt consistency checks. Not a long-term growth exponent.',
     'definition':'K_delta=sum_i m_i |v_i-v_circular(x_i)|^2/2, fixed initial material regions.',
     'script_sha256':reader.sha(Path(__file__)),'regions':{},'accepted_step_sources':[]}
for name,mask in masks.items():
    p=power[name];h=hg[name]['work_on_actual_velocity_departure']
    assert abs(energy0[mask].sum()-p['kinetic_energy_relative_to_circular'])<1e-16
    out['regions'][name]={'cells':int(mask.sum()),'pressure':p['pressure_work'],'TCF':p['TCF_work']['saved_chi'],
                         'HG':h,'differential_rotation':p['differential_rotation_work'],
                         'full_instantaneous_power':p['total_without_HG']['saved_chi']+h,
                         'chi_zero_instantaneous_power':p['total_without_HG']['chi_zero']+h,
                         'finite_dt_checks':[]}
for case in ['t6_full','t6_half','t6_quarter']:
    d=base/'Codex_early_growth_budget_20261017'/case
    assert (d/'exit.code').read_text().strip()=='0'
    assert reader.sha(d/'khout.000656.dat')==reader.sha(root/'khout.000656.dat')
    a,t=reader.read(d/'khout.000657.dat');assert np.array_equal(a['mass'],start['mass'])
    rate=(cell_energy(a)-energy0)/(t-t0)
    out['accepted_step_sources'].append({'run':str(d),'dt':t-t0,'end_sha256':reader.sha(d/'khout.000657.dat')})
    for name,mask in masks.items():
        exact=out['regions'][name]['full_instantaneous_power'];measured=float(rate[mask].sum())
        out['regions'][name]['finite_dt_checks'].append({'dt':t-t0,'measured_mean_power':measured,
                                                       'difference_from_instantaneous':measured-exact})
for name,row in out['regions'].items():
    z=[v['measured_mean_power'] for v in row['finite_dt_checks']]
    row['extrapolated_power_coarse']=2*z[1]-z[0]
    row['extrapolated_power_fine']=2*z[2]-z[1]
    row['extrapolated_error_coarse']=row['extrapolated_power_coarse']-row['full_instantaneous_power']
    row['extrapolated_error_fine']=row['extrapolated_power_fine']-row['full_instantaneous_power']
with (root/'full_power.json').open('x') as f:json.dump(out,f,indent=2);f.write('\n')
print(json.dumps(out,indent=2))

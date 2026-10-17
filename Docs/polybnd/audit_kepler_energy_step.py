"""Active gas-plus-potential energy balance for a saved split Kepler trial.

Accounts for gas boundary ledger and independently integrated potential mass
exchange with the frozen reservoir. No energy correction or time replay.
"""
import argparse,hashlib,json,os
from pathlib import Path
import numpy as np
import moment_ale_evolve as base
from physical_field_accuracy_audit import environment
from check_kepler_moment_preflight import install
from adaptive_coupled_mass import restore
from gravity_coupled_2d import hydro_dense_state
from smooth_mass_path import evaluate_smooth_path
from kepler_disk_reference import polar_rule,frozen,CENTER
from kepler_spline_gravity import potential
from check_kepler_dense_step import atomic,simpson


def phi(x):return potential(x,center=CENTER,rc=.5)


def configure(cfg):
    os.environ.update(environment(cfg));os.environ.update(LAGFORCE_VECTOR_KINETIC_ROOT='1',
        LAGFORCE_KINETIC_SMOOTHING='1',LAGFORCE_KINETIC_SMOOTHING_LENGTH='1',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star');install()


def initial_key(source,configuration,order):
    return dict(source_sha256=hashlib.sha256(source).hexdigest(),configuration_sha256=hashlib.sha256(configuration).hexdigest(),
        source_manifest_sha256=hashlib.sha256(Path(__file__).with_name('SHA256SUMS').read_bytes()).hexdigest(),order=order)


def balance(gas_initial,gas_final,potential_initial,potential_final,gas_exchange,potential_exchange):
    return (gas_final-gas_initial)+(potential_final-potential_initial)-gas_exchange-potential_exchange


def active_potential(cp,mass,mask,order,path,identity):
    a=np.asarray(cp['state'])[:-4].reshape(-1,20)
    cells,_=base.geometry(a[:,:2],a[:,2],cp['t'])
    saved=dict(identity=identity,next_cell=0,energy=0.,maximum_mass_quadrature_error=0.)
    if Path(path).exists():
        saved=json.loads(Path(path).read_text())
        if saved['identity']!=identity:raise ValueError('potential cache identity changed')
    for i in range(saved['next_cell'],len(cells)):
        if not mask[i]:
            density=base.decode(cells[i],np.r_[mass[i],a[i,3:]])['density']
            x,w=polar_rule(cells[i],order,radii=(.25,.5));dm=w*density.value(x)
            saved['energy']+=float(dm@phi(x))
            saved['maximum_mass_quadrature_error']=max(saved['maximum_mass_quadrature_error'],float(abs(dm.sum()/mass[i]-1)))
        if (i+1)%32==0 or i+1==len(cells):
            saved['next_cell']=i+1;atomic(path,saved)
            print('potential quadrature',order,i+1,'/',len(cells),flush=True)
    return saved


def prepare_initial(seed,folder):
    root=Path(seed);folder=Path(folder);folder.mkdir(exist_ok=True)
    source=(root/'source.json').read_bytes();configuration=(root/'configuration.json').read_bytes()
    cp=json.loads(source);configure(json.loads(configuration));mass=np.asarray(cp['mass'])
    mask=frozen(cp['t'],np.asarray(cp['state'])[:-4].reshape(-1,20)[:,:2])
    for order in (16,24):
        key=initial_key(source,configuration,order)
        result=active_potential(cp,mass,mask,order,folder/f'initial_{order}_progress.json',key)
        atomic(folder/f'initial_{order}_result.json',result)


def run(seed,initial_cache=None):
    root=Path(seed);names=('source.json','trial.json','configuration.json','nodes.json')
    data={n:(root/n).read_bytes() for n in names}
    identity=dict(inputs={n:hashlib.sha256(v).hexdigest() for n,v in data.items()},
        source_manifest_sha256=hashlib.sha256(Path(__file__).with_name('SHA256SUMS').read_bytes()).hexdigest())
    if Path('identity.json').exists():
        if json.loads(Path('identity.json').read_text())!=identity:raise ValueError('energy audit identity changed')
    else:atomic('identity.json',identity)
    cp=json.loads(data['source.json']);trial=json.loads(data['trial.json']);end=trial['endpoint']
    cfg=json.loads(data['configuration.json']);record=restore(trial['record']);mass=np.asarray(cp['mass'])
    configure(cfg)
    np.testing.assert_array_equal(mass,end['mass'])
    mask=np.asarray(record['frozen_first'],dtype=bool)
    for state in (cp,end):
        sites=np.asarray(state['state'])[:-4].reshape(-1,20)[:,:2]
        np.testing.assert_array_equal(frozen(state['t'],sites),mask)
    samples=json.loads(Path('exchange_nodes.json').read_text()) if Path('exchange_nodes.json').exists() else []
    reference=json.loads(data['nodes.json'])
    if len(reference)!=17:raise ValueError('requires completed 17-node seed audit')
    for node in range(len(samples),17):
        s=node/16;y,L,v=hydro_dense_state(record,mass,s);a=v[:-4].reshape(-1,20)
        np.testing.assert_array_equal(frozen(record['t']+s*record['dt'],y[:-4].reshape(-1,20)[:,:2]),mask)
        flux,d=evaluate_smooth_path(record['t']+s*record['dt'],y,L,mass,(a[:,:2],a[:,2]),
            frozen_mask=mask,potential=phi)
        np.testing.assert_allclose(flux,reference[node]['flux'],atol=2e-13,rtol=1e-10)
        samples.append(d['active_potential_mass_exchange']);atomic('exchange_nodes.json',samples)
        print('potential exchange node',node,'/16',flush=True)
    exchange=float(simpson(np.asarray(samples),record['dt']))
    temporal_difference=float(abs(exchange-simpson(np.asarray(samples)[::2],record['dt'])))
    old=np.asarray(cp['state']);new=np.asarray(end['state']);a=old[:-4].reshape(-1,20);b=new[:-4].reshape(-1,20)
    np.testing.assert_array_equal(a[mask],b[mask])
    gas0=float(a[~mask,19].sum());gas1=float(b[~mask,19].sum());gas_exchange=float(new[-2]-old[-2])
    thermal=float((a[~mask,19]-.5*(a[~mask,15]+a[~mask,18])).sum())
    rows=[]
    for order in (16,24):
        energies=[]
        for label,state in [('initial',cp),('final',end)]:
            if label=='initial' and initial_cache is not None:
                cached=json.loads((Path(initial_cache)/f'initial_{order}_result.json').read_text())
                if cached['identity']!=initial_key(data['source.json'],data['configuration.json'],order):raise ValueError('shared initial potential identity mismatch')
                energies.append(cached)
            else:energies.append(active_potential(state,mass,mask,order,f'{label}_{order}_progress.json',dict(identity=identity,state=label,order=order)))
        residual=balance(gas0,gas1,energies[0]['energy'],energies[1]['energy'],gas_exchange,exchange)
        row=dict(order=order,initial_potential=energies[0]['energy'],final_potential=energies[1]['energy'],
            residual=residual,relative_to_initial_thermal=abs(residual)/thermal,
            maximum_mass_quadrature_error=max(x['maximum_mass_quadrature_error'] for x in energies))
        rows.append(row);atomic(f'result_{order}.json',row)
    report=dict(t=cp['t'],dt=record['dt'],gas_energy_change=gas1-gas0,gas_boundary_exchange=gas_exchange,
        potential_boundary_exchange=exchange,potential_exchange_nested_difference=temporal_difference,
        gas_ledger_gravity_identity_error=abs(gas1-gas0-gas_exchange-trial['diagnostics']['gravity_work']),
        initial_active_thermal=thermal,rows=rows,potential_change_quadrature_difference=abs(rows[1]['residual']-rows[0]['residual']),full_validation=False)
    atomic('report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('seed');p.add_argument('--prepare-initial');p.add_argument('--initial-cache');a=p.parse_args()
    if a.prepare_initial:prepare_initial(a.seed,a.prepare_initial)
    else:run(a.seed,a.initial_cache)

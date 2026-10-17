"""Semi-discrete potential/work compatibility at a saved hydro stage.

Symmetric directional differences of reconstructed potential energy use the
actual saved hydro derivative. No new integration step or energy correction.
"""
import argparse, hashlib, json, os
from pathlib import Path
import numpy as np
import moment_ale_evolve as base
from adaptive_coupled_mass import restore
from gravity_coupled_2d import hydro_dense_state
from audit_kepler_energy_step import configure, active_potential, phi, atomic
from kepler_disk_reference import CENTER, polar_rule
from kepler_spline_gravity import acceleration
from kinetic_factor_rhs import ALPHA
from quadratic_velocity_reconstruction import prepare, velocity_samples
from smooth_mass_path import evaluate_smooth_path


def run(seed):
    root=Path(seed)
    files={n:(root/n).read_bytes() for n in ('source.json','configuration.json','trial.json')}
    identity={n:hashlib.sha256(v).hexdigest() for n,v in files.items()}
    identity['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    cp=json.loads(files['source.json']);configure(json.loads(files['configuration.json']))
    record=restore(json.loads(files['trial.json'])['record']);mass=np.asarray(cp['mass'])
    y,L,dy=hydro_dense_state(record,mass,0.)
    mask=np.asarray(record['frozen_first'],dtype=bool);a=y[:-4].reshape(-1,20)
    cells,edges=base.geometry(a[:,:2],a[:,2],record['t'])
    states=[base.decode(p,np.r_[m,row[3:]]) for p,m,row in zip(cells,mass,a)]
    if os.environ.get('LAGFORCE_QUADRATIC_VELOCITY')=='1':prepare(states,cells,edges,mass,L)
    work_native=0.;work_polar=0.
    for i,s in enumerate(states):
        if mask[i]:continue
        density=s['density']
        for mode in ('native','polar'):
            if mode=='native':x=density.points;dm=density.weighted_mass
            else:
                x,w=polar_rule(cells[i],16,radii=(.25,.5));dm=w*density.value(x)
            v,_=velocity_samples(s,L[i],mass[i],x,ALPHA)
            work=float(dm@np.sum(v.mean(1)*acceleration(x,center=CENTER,rc=.5),axis=1))
            if mode=='native':work_native+=work
            else:work_polar+=work
    rates=dy[:-4].reshape(-1,20)
    _,diag=evaluate_smooth_path(record['t'],y,L,mass,(rates[:,:2],rates[:,2]),frozen_mask=mask,potential=phi)
    boundary=diag['active_potential_mass_exchange'];rows=[]
    for h in (1e-5,5e-6):
        values=[]
        for sign in (-1,1):
            label=f'h{h}_sign{sign}'
            perturbed=dict(t=record['t']+sign*h,state=y+sign*h*dy)
            result=active_potential(perturbed,mass,mask,16,label+'.json',dict(identity=identity,h=h,sign=sign))
            values.append(result['energy'])
        slope=(values[1]-values[0])/(2*h)
        rows.append(dict(h=h,potential_directional_derivative=slope,
            residual_native=slope+work_native-boundary,
            residual_polar=slope+work_polar-boundary))
    report=dict(identity=identity,t=record['t'],gravity_work_rate_native=work_native,
        gravity_work_rate_polar=work_polar,potential_boundary_rate=boundary,rows=rows,
        full_validation=False,scope='single-stage directional consistency, fixed original masses')
    atomic('report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('seed');run(p.parse_args().seed)

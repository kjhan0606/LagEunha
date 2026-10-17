"""Check all degree <=2 potential identities on an actual saved hydro stage.

Potential moments and gravity work use stored q,Q,p,D exactly. Boundary
potential flux is independently evaluated with the actual shared mass flux.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import restore
from audit_kepler_energy_step import configure,atomic
from gravity_coupled_2d import hydro_dense_state
from smooth_mass_path import evaluate_smooth_path


def potential_basis(x):
    x=np.asarray(x)
    return np.column_stack([np.ones(len(x)),x[:,0],x[:,1],.5*x[:,0]**2,x[:,0]*x[:,1],.5*x[:,1]**2])


def moment_rates(a,da):
    """dU/dt and -integral rho v.grad(Phi) for six fixed potentials."""
    udot=np.column_stack([np.zeros(len(a)),da[:,3],da[:,4],.5*da[:,5],.5*(da[:,6]+da[:,7]),.5*da[:,8]])
    work=np.column_stack([np.zeros(len(a)),-a[:,9],-a[:,10],-a[:,11],-a[:,12]-a[:,13],-a[:,14]])
    return udot,work


def run(seed):
    root=Path(seed);names=('source.json','trial.json','configuration.json')
    data={n:(root/n).read_bytes() for n in names}
    configure(json.loads(data['configuration.json']))
    mass=np.asarray(json.loads(data['source.json'])['mass'])
    record=restore(json.loads(data['trial.json'])['record'])
    y,L,dy=hydro_dense_state(record,mass,0.)
    a=y[:-4].reshape(-1,20);da=dy[:-4].reshape(-1,20)
    mask=np.asarray(record['frozen_first'],dtype=bool)
    udot,work=moment_rates(a,da);udot=udot[~mask].sum(0);work=work[~mask].sum(0)
    rows=[]
    for k,name in enumerate(['1','x','y','x2/2','xy','y2/2']):
        flux,d=evaluate_smooth_path(record['t'],y,L,mass,(da[:,:2],da[:,2]),
            frozen_mask=mask,potential=lambda x:potential_basis(x)[:,k])
        exchange=d['active_potential_mass_exchange']
        rows.append(dict(potential=name,potential_rate=float(udot[k]),gravity_work_rate=float(work[k]),
            boundary_rate=exchange,residual=float(udot[k]+work[k]-exchange),
            maximum_active_mass_rate=float(np.max(abs(flux[~mask])))))
    report=dict(identity={n:hashlib.sha256(v).hexdigest() for n,v in data.items()},
        rows=rows,maximum_residual=max(abs(r['residual']) for r in rows),full_validation=False)
    atomic('report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('seed');run(p.parse_args().seed)

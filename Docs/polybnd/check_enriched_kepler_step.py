"""Independent saved-state trial of coupled potential-moment transport."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import moment_ale_evolve as base
from audit_kepler_energy_step import configure,phi,atomic
from kepler_disk_reference import polar_rule,frozen,CENTER
from kepler_spline_gravity import acceleration
from potential_moment_coupled_2d import initial_potential,step,provider


def run(seed,dt,gravity_scheme='finite'):
    root=Path(seed);data={n:(root/n).read_bytes() for n in ('source.json','configuration.json')}
    identity={n:hashlib.sha256(v).hexdigest() for n,v in data.items()};identity['dt']=dt
    identity['gravity_scheme']=gravity_scheme
    identity['source_manifest_sha256']=hashlib.sha256(Path(__file__).with_name('SHA256SUMS').read_bytes()).hexdigest()
    if Path('identity.json').exists():
        if json.loads(Path('identity.json').read_text())!=identity:raise ValueError('trial identity changed')
    else:atomic('identity.json',identity)
    cp=json.loads(data['source.json']);configure(json.loads(data['configuration.json']))
    y=np.asarray(cp['state']);mass=np.asarray(cp['mass']);L=[np.asarray(a) for a in cp['factors']];t=cp['t']
    quad=lambda p:polar_rule(p,16,radii=(.25,.5));g=lambda x:acceleration(x,center=CENTER,rc=.5)
    if not Path('initial_potential.json').exists():atomic('initial_potential.json',initial_potential(t,y,mass,phi,quad))
    U=np.asarray(json.loads(Path('initial_potential.json').read_text()))
    print('initial potential ready',flush=True)
    if not Path('trial.json').exists():
        out,F,newU,d=step(t,y,mass,L,U,dt,phi,g,quad,frozen,return_record=True,gravity_scheme=gravity_scheme)
        record=d.pop('hydro_record')
        atomic('trial.json',dict(t=t+dt,state=out,factors=F,potential_moments=newU,diagnostics=d,mass=mass,hydro_record=record))
    trial=json.loads(Path('trial.json').read_text());out=np.asarray(trial['state']);newU=np.asarray(trial['potential_moments']);d=trial['diagnostics']
    a=y[:-4].reshape(-1,20);b=out[:-4].reshape(-1,20);mask=frozen(t,a[:,:2])
    np.testing.assert_array_equal(mass,trial['mass']);np.testing.assert_array_equal(a[mask],b[mask]);np.testing.assert_array_equal(U[mask],newU[mask])
    for i in np.flatnonzero(mask):np.testing.assert_array_equal(L[i],trial['factors'][i])
    gas=float((b[~mask,19]-a[~mask,19]).sum());du=float((newU[~mask]-U[~mask]).sum());gas_boundary=float(out[-2]-y[-2])
    residual=gas+du-gas_boundary-d['potential_boundary_exchange']
    cells,_=base.geometry(b[:,:2],b[:,2],t+dt);get=provider(newU,phi,quad,mask);mass_error=0.;potential_error=0.
    for i,p in enumerate(cells):
        if mask[i]:continue
        density=get(i,p,np.r_[mass[i],b[i,3:]]);x,w=polar_rule(p,24,radii=(.25,.5));dm=w*density.value(x)
        mass_error=max(mass_error,float(abs(dm.sum()/mass[i]-1)));potential_error=max(potential_error,float(abs(dm@phi(x)-newU[i])))
    report=dict(t=t,dt=dt,diagnostics=d,energy_residual=residual,gas_energy_change=gas,potential_change=du,
        gas_boundary_exchange=gas_boundary,independent_mass_error=mass_error,independent_potential_moment_error=potential_error,
        full_validation=False,endpoint_accepted=False,scope='single coupled saved-state trial with extra potential-density moment')
    atomic('report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('seed');p.add_argument('--dt',type=float,required=True);p.add_argument('--gravity-scheme',choices=('finite','midpoint'),default='finite');a=p.parse_args();run(a.seed,a.dt,a.gravity_scheme)

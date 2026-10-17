"""Small-amplitude 3D standing acoustic mode; linear-reference O(A^2) limit.

This first physical evolution gate uses explicitly linear pressure and affine
velocity. Fixed primary masses and stage constraints are checked; continuous
moving-face mass-path certification is still required separately.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from power_cells_3d import Geometry3D
from moment_density_3d import quadrature
from moment_state_3d import pack,decode
from moment_evolution_3d import step

GAMMA=5/3;AMPLITUDE=1e-5;C=np.sqrt(GAMMA)


def reference(x,t):
    mode=np.cos(2*np.pi*x[:,0])*np.cos(2*np.pi*C*t)
    rho=1+AMPLITUDE*mode;P=1+GAMMA*AMPLITUDE*mode;v=np.zeros_like(x)
    v[:,0]=AMPLITUDE*C*np.sin(2*np.pi*x[:,0])*np.sin(2*np.pi*C*t)
    return rho,P,v


def score(geo,sites,weights,raw,t,pressure_model="linear"):
    polys,_,_=geo.build(sites,weights);states=[decode(poly,r) for poly,r in zip(polys,raw)]
    centers=np.array([s['density'].origin for s in states]);_,faces,_=geo.build(sites,weights);neighbors=[set() for _ in sites]
    for f in faces:
        if f['j']>=0:neighbors[f['i']].add(f['j']);neighbors[f['j']].add(f['i'])
    pressure_polynomials=None
    if pressure_model=="quadratic":
        from pressure_reconstruction_3d import prepare,value
        pressure_polynomials,_=prepare(states,polys)
    error=np.zeros(3);minimum=np.full(2,np.inf)
    for i,(poly,s) in enumerate(zip(polys,states)):
        ids=list(neighbors[i]);dr=centers[ids]-centers[i];d=np.linalg.norm(dr,axis=1);dp=np.array([states[j]['P']-s['P'] for j in ids]);g=np.linalg.lstsq(dr/d[:,None],dp/d,rcond=None)[0]
        x,w=quadrature(poly,8);rho,P,v=reference(x,t);rn=s['density'].value(x);pn=s['P']+(x-centers[i])@g;vn=s['v']+(x-s['c'])@s['B'].T
        if pressure_polynomials is not None:pn=value(pressure_polynomials[i],x)
        error+=np.array([w@abs(rn-rho)/AMPLITUDE,w@np.linalg.norm(vn-v,axis=1)/(AMPLITUDE*C),w@abs(pn-P)/(AMPLITUDE*GAMMA)])
        minimum=np.minimum(minimum,[rn.min(),pn.min()])
    return dict(density_L1=float(error[0]),velocity_L1=float(error[1]),pressure_L1=float(error[2]),min_density=float(minimum[0]),min_pressure=float(minimum[1]))


def load_resume(n,dt,end,pressure_model):
    cp=json.loads(Path('checkpoint.json').read_text());ic=json.loads(Path('initial_state.json').read_text())
    for key,value in dict(n=n,dt=dt,end=end,pressure_model=pressure_model).items():
        if cp[key]!=value or ic[key]!=value:raise ValueError('resume configuration mismatch: '+key)
    if not 0.<cp['t']<=end:raise ValueError('invalid checkpoint time')
    mass=np.asarray(cp['mass']);raw=np.asarray(cp['raw'])
    if not np.array_equal(mass,np.asarray(ic['mass'])) or not np.array_equal(raw[:,0],mass):raise ValueError('resume changed original masses')
    for key in ('sites','weights','raw','factors','ledger'):
        if not np.isfinite(np.asarray(cp[key])).all():raise ValueError('nonfinite checkpoint: '+key)
    return cp,ic


def run(library,n,dt,end,pressure_model="linear",resume=False):
    geo=Geometry3D(library);resumed=False;peak_is_bound=False
    if resume:
        cp,ic=load_resume(n,dt,end,pressure_model)
        sites=np.asarray(cp['sites']);weights=np.asarray(cp['weights']);raw=np.asarray(cp['raw']);L=np.asarray(cp['factors'])
        mass=np.asarray(cp['mass']);original=np.asarray(ic['raw']);ledger=np.asarray(cp['ledger']);t=cp['t']
        initial=json.loads(Path('initial_score.json').read_text());resumed=True
        # Legacy checkpoints omitted the measured peak. Each accepted stage
        # enforced this bound; report it as a bound, never a measured maximum.
        peak=cp.get('stage_mass_flux_peak',2e-13);peak_is_bound=cp.get('stage_mass_flux_is_bound','stage_mass_flux_peak' not in cp)
        start=time.monotonic()
    else:
        if Path('checkpoint.json').exists() or Path('initial_state.json').exists():raise RuntimeError('existing acoustic state: use explicit resume')
        rng=np.random.default_rng(1021+n);grid=np.meshgrid(*([(np.arange(n)+.5)/n]*3),indexing='ij');sites=np.column_stack([a.ravel() for a in grid])+rng.uniform(-.05,.05,(n**3,3))/n;weights=np.zeros(len(sites));polys,_,_=geo.build(sites,weights);raw=[];L=np.zeros((len(sites),3,3))
        for poly in polys:
            x,w=quadrature(poly,8);rho,P,v=reference(x,0.);dm=w*rho
            raw.append(pack(dm.sum(),dm@x,x.T@(dm[:,None]*x),np.zeros(3),np.zeros((3,3)),np.zeros((3,3)),w@P/(GAMMA-1)))
        raw=np.array(raw);mass=raw[:,0].copy();original=raw.copy();initial=score(geo,sites,weights,raw,0.,pressure_model);ledger=np.zeros(7);peak=0.;t=0.;start=time.monotonic()
        Path('initial_state.json').write_text(json.dumps(dict(t=0.,sites=sites.tolist(),weights=weights.tolist(),raw=raw.tolist(),factors=L.tolist(),mass=mass.tolist(),n=n,dt=dt,end=end,pressure_model=pressure_model)))
        Path('initial_score.json').write_text(json.dumps(initial,indent=2))
    bounds=lambda t:np.array([[0.,0.,0.],[1.,1.,1.]])
    while t<end-1e-15:
        h=min(dt,end-t);sites,weights,raw,L,d=step(geo,t,sites,weights,raw,L,h,bounds,lambda p:np.zeros(3),pressure_reconstruction=pressure_model);t+=h;ledger+=d['ledger_increment'];peak=max(peak,d['stage_mass_flux']);assert np.array_equal(mass,raw[:,0])
        cp=dict(stage_mass_flux_peak=peak,stage_mass_flux_is_bound=peak_is_bound,t=t,sites=sites.tolist(),weights=weights.tolist(),raw=raw.tolist(),factors=L.tolist(),ledger=ledger.tolist(),mass=mass.tolist(),n=n,dt=dt,end=end,pressure_model=pressure_model)
        Path('checkpoint.tmp').write_text(json.dumps(cp));Path('checkpoint.tmp').replace('checkpoint.json');print('3D acoustic',n,t,flush=True)
    final=score(geo,sites,weights,raw,t,pressure_model);result=dict(resumed=resumed,stage_mass_flux_is_bound=peak_is_bound,pressure_model=pressure_model,completed=True,n=n,t=t,dt=dt,amplitude=AMPLITUDE,initial=initial,final=final,stage_mass_flux=peak,mass_exact=True,energy_ledger_error=float(abs(raw[:,34].sum()-original[:,34].sum()-ledger[3])),physical_goal=.05,physical_passed=max(final[k] for k in ('density_L1','velocity_L1','pressure_L1'))<.05,dense_path_certified=False,seconds=time.monotonic()-start)
    Path('acoustic_3d_report.json').write_text(json.dumps(result,indent=2));return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('library');p.add_argument('--n',type=int,default=3);p.add_argument('--dt',type=float,default=.002);p.add_argument('--end',type=float,default=.02);p.add_argument('--pressure',choices=['linear','quadratic'],default='linear');p.add_argument('--resume',action='store_true');a=p.parse_args();print(json.dumps(run(a.library,a.n,a.dt,a.end,a.pressure,a.resume)),flush=True)

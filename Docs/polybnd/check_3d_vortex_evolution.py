"""Nonlinear stationary Euler vortex extruded in z on actual 3D moving cells.

No gravity/steering, original fixed masses, affine velocity and quadratic
pressure. Save every actual dense path losslessly for subsequent mass audits.
Reference fields are used only in initialization/scoring.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from steady_vortex_exact import fields as fields2d
from power_cells_3d import Geometry3D
from moment_density_3d import quadrature
from moment_gravity_kick import moments
from moment_state_3d import pack,decode
from moment_evolution_3d import step,rebuild
from pressure_reconstruction_3d import prepare,value
from dense_record_3d import save,load


def atomic(path,data):
    p=Path(path);tmp=Path(str(p)+'.tmp');tmp.write_text(json.dumps(data));tmp.replace(p)


def fields(x,amplitude):
    rho,v,P=fields2d(x[:,:2],amplitude=amplitude)
    return rho,np.column_stack([v,np.zeros(len(v))]),P


def initial_raw(poly,amplitude,order):
    x,w=quadrature(poly,order);rho,v,P=fields(x,amplitude);d=moments(x,w*rho,v[:,None,:],w@P/(5/3-1))
    return pack(*(d[k] for k in ('m','q','Q','p','D','Pi','E'))),d['residual_factor']


def score(geo,sites,weights,raw,amplitude):
    polys,_,_=geo.build(sites,weights);states=[decode(p,r) for p,r in zip(polys,raw)];pressure,_=prepare(states,polys);rows=[]
    for order in (8,12):
        errors=np.zeros(3);denom=0.;minimum=np.full(2,np.inf);max_vz=0.
        for i,(poly,s) in enumerate(zip(polys,states)):
            x,w=quadrature(poly,order);rho,v,P=fields(x,amplitude);rn=s['density'].value(x);vn=s['v']+(x-s['c'])@s['B'].T;pn=value(pressure[i],x)
            errors+=w@np.column_stack([abs(rn-rho),np.linalg.norm(vn-v,axis=1)/amplitude,abs(pn-P)/amplitude**2]);denom+=w@np.linalg.norm(v,axis=1)
            minimum=np.minimum(minimum,[rn.min(),pn.min()]);max_vz=max(max_vz,float(np.max(abs(vn[:,2]))))
        rows.append(dict(order=order,density_L1=float(errors[0]),velocity_L1=float(errors[1]),pressure_L1=float(errors[2]),velocity_reference_relative_L1=float(errors[1]*amplitude/denom),min_density=float(minimum[0]),min_pressure=float(minimum[1]),max_vertical_velocity_sampled=max_vz))
    return rows


def run(library,n,dt,end,amplitude):
    cfg=dict(n=n,dt=dt,end=end,amplitude=amplitude,pressure_model='quadratic',velocity_model='affine')
    if Path('configuration.json').exists():
        if json.loads(Path('configuration.json').read_text())!=cfg:raise ValueError('configuration mismatch')
    else:atomic('configuration.json',cfg)
    geo=Geometry3D(library);start=time.monotonic()
    if not Path('initial_state.json').exists():
        rng=np.random.default_rng(1021+n);grid=np.meshgrid(*([(np.arange(n)+.5)/n]*3),indexing='ij');sites=np.column_stack([a.ravel() for a in grid])+rng.uniform(-.05,.05,(n**3,3))/n;weights=np.zeros(len(sites));polys,_,_=geo.build(sites,weights);raw=[];factors=[];reference_error=0.
        for i,p in enumerate(polys):
            low,_=initial_raw(p,amplitude,10);r,L=initial_raw(p,amplitude,14);error=float(np.max(abs(low-r))/np.max(abs(r)));reference_error=max(reference_error,error)
            if error>1e-5:raise RuntimeError('unresolved 3D vortex IC quadrature')
            raw.append(r);factors.append(L)
        raw=np.asarray(raw);factors=np.asarray(factors);original=raw.copy();raw=rebuild(raw,factors)
        ic=dict(t=0.,sites=sites.tolist(),weights=weights.tolist(),raw=raw.tolist(),factors=factors.tolist(),mass=raw[:,0].tolist(),ledger=[0.]*7,steps=0,stage_mass_peak=0.,configuration=cfg,reference_quadrature_difference=reference_error,Pi_initial_completion=float(np.max(abs(raw-original))))
        atomic('initial_state.json',ic)
    ic=json.loads(Path('initial_state.json').read_text())
    if not Path('initial_score.json').exists():atomic('initial_score.json',score(geo,np.asarray(ic['sites']),np.asarray(ic['weights']),np.asarray(ic['raw']),amplitude))
    cp=json.loads(Path('checkpoint.json').read_text()) if Path('checkpoint.json').exists() else ic
    if cp['configuration']!=cfg:raise ValueError('checkpoint configuration mismatch')
    mass=np.asarray(ic['mass']);sites=np.asarray(cp['sites']);weights=np.asarray(cp['weights']);raw=np.asarray(cp['raw']);L=np.asarray(cp['factors']);ledger=np.asarray(cp['ledger']);t=cp['t'];count=cp['steps'];peak=cp['stage_mass_peak'];np.testing.assert_array_equal(raw[:,0],mass)
    bounds=lambda t:np.array([[0.,0.,0.],[1.,1.,1.]])
    while t<end-1e-15:
        h=min(dt,end-t);path=Path(f'step_{count+1:05d}.npz')
        if path.exists():
            record,endpoint=load(path)
            if record['t']!=t or record['dt']!=h:raise ValueError('pending dense record time mismatch')
            for k,v in dict(sites=sites,weights=weights,raw=raw,factors=L).items():np.testing.assert_array_equal(record[k],v)
            next_sites=endpoint['sites'];next_weights=endpoint['weights'];next_raw=endpoint['raw'];nextL=endpoint['factors'];increment=endpoint['ledger'];stage_mass=float(endpoint['stage_mass'])
        else:
            next_sites,next_weights,next_raw,nextL,d=step(geo,t,sites,weights,raw,L,h,bounds,lambda p:np.zeros(3),pressure_reconstruction='quadratic',return_record=True)
            increment=d['ledger_increment'];stage_mass=d['stage_mass_flux'];save(path,d['dense_record'],dict(sites=next_sites,weights=next_weights,raw=next_raw,factors=nextL,ledger=increment,stage_mass=stage_mass))
        np.testing.assert_array_equal(next_raw[:,0],mass)
        if stage_mass>=2e-13:raise RuntimeError('3D stage mass gate failed')
        sites,weights,raw,L=next_sites,next_weights,next_raw,nextL;t+=h;count+=1;ledger+=increment;peak=max(peak,stage_mass)
        cp=dict(t=t,sites=sites.tolist(),weights=weights.tolist(),raw=raw.tolist(),factors=L.tolist(),mass=mass.tolist(),ledger=ledger.tolist(),steps=count,stage_mass_peak=peak,configuration=cfg);atomic('checkpoint.json',cp);print('3D nonlinear vortex',n,t,flush=True)
    final=score(geo,sites,weights,raw,amplitude);passed=all(max(r[k] for k in ('density_L1','velocity_L1','pressure_L1'))<.05 for r in final)
    report=dict(completed=True,n=n,t=t,dt=dt,amplitude=amplitude,initial=json.loads(Path('initial_score.json').read_text()),final=final,physical_goal=.05,physical_passed=passed,mass_exact=True,stage_mass_peak=peak,energy_ledger_error=float(abs(raw[:,34].sum()-np.asarray(ic['raw'])[:,34].sum()-ledger[3])),dense_records_saved=count,dense_path_certified=False,seconds=time.monotonic()-start)
    atomic('vortex_3d_report.json',report);print(json.dumps(report),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('library');p.add_argument('--n',type=int,required=True);p.add_argument('--dt',type=float,default=.002);p.add_argument('--end',type=float,default=.02);p.add_argument('--amplitude',type=float,default=.3);a=p.parse_args();run(a.library,a.n,a.dt,a.end,a.amplitude)

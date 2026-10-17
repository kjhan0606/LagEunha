"""Independent stationary Euler Reynolds-rate comparison, on moving cells."""
import json,time,hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from adaptive_coupled_mass import atomic
from steady_vortex_exact import physical_raw,fields,parameters,GAMMA
from steady_vortex_evolution import initial,score
from simple_wave_evolution import install_domain
from simple_wave_rhs_consistency import exact_rates
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


def time_fields(points,t,amplitude,with_time_derivatives=False):
    rho,v,P=fields(points,t,amplitude)
    if with_time_derivatives:return rho,v,P,np.zeros_like(rho),np.zeros_like(v),np.zeros_like(P)
    return rho,v,P


def inspect(entry):
    started=time.monotonic();install_domain();amp=.5;_,drop=parameters(amp)
    n=entry['n'] if isinstance(entry,dict) else entry
    if isinstance(entry,dict):
        from saved_curved_contact import checked_json
        source=checked_json(entry['initial'],entry['initial_sha256'])
        if source['physical_initial_condition']!=dict(kind='steady_vortex',amplitude=amp,n=n):raise ValueError('vortex source identity mismatch')
        mass=np.asarray(source['mass']);state=np.asarray(source['state']);L=[np.asarray(f) for f in source['factors']]
        adjustment=source.get('roundoff_Pi_completion',0.)
    else:mass,state,L,adjustment=initial(n,amp)
    dy,diag=kinetic.rhs(0.,state,mass,covariance_factors=L,face_order=4)
    atomic('n'+str(n)+'_initial_state.json',dict(n=n,t=0.,mass=mass,state=state,factors=L,
        physical_initial_condition=dict(kind='steady_vortex',amplitude=amp,n=n),roundoff_Pi_completion=adjustment))
    atomic('n'+str(n)+'_source_rhs.json',dict(t=0.,state_sha256=hashlib.sha256(state.tobytes()).hexdigest(),
        mass_sha256=hashlib.sha256(mass.tobytes()).hexdigest(),dy=dy,diag=diag))
    a=state[:-4].reshape(-1,20);d=dy[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],0.)
    def reference(order):return exact_rates(cells,edges,a[:,:2],d[:,:2],d[:,2],0.,amp,order,time_fields,physical_raw)
    r,area,area_rate=reference(16);fine,_,_=reference(24)
    P=(GAMMA-1)*(a[:,19]-.5*(a[:,15]+a[:,18]))/area
    pn=(GAMMA-1)*(d[:,19]-.5*(d[:,15]+d[:,18]))/area-P*area_rate/area
    pr=(GAMMA-1)*(r[:,17]-.5*(r[:,13]+r[:,16]))/area-P*area_rate/area
    vr=(r[:,7:9]-a[:,9:11]/mass[:,None]*r[:,0,None])/mass[:,None]
    errors=dict(pressure_mean_rate=float(area@abs(pn-pr)/drop),velocity_mean_rate=float(area@np.linalg.norm(d[:,9:11]/mass[:,None]-vr,axis=1)/amp),
        density_mass_rate=float(abs(r[:,0]).sum()),momentum_rate=float(np.linalg.norm(d[:,9:11]-r[:,7:9],axis=1).sum()/amp),
        energy_rate=float(abs(d[:,19]-r[:,17]).sum()/drop))
    result=dict(n=n,errors=errors,initial_fields=score(0.,state,mass,amp,L),reference_rate_quadrature_difference=float(np.max(np.sum(abs(r-fine),axis=0))),
        source_mass_flux=diag['mass_flux'],velocity_reconstruction=diag['velocity_reconstruction'],pressure_reconstruction=diag.get('pressure_reconstruction','linear'),velocity_budget_diagnostics={k:v for k,v in diag.items() if k.startswith(('velocity_quadratic','pressure_quadratic'))},seconds=time.monotonic()-started,
        global_exact_rates=dict(mass=float(r[:,0].sum()),momentum=r[:,7:9].sum(0),energy=float(r[:,17].sum()),angular_momentum=float((r[:,10]-r[:,11]).sum())),
        note='Stationary Euler reference has zero Eulerian derivatives, but actual moving-cell reference rates include Reynolds boundary motion. Exact field never enters RHS.')
    atomic('n'+str(n)+'_vortex_rhs.json',result);print(result,flush=True);return result


def run():
    config=Path('vortex_rhs_config.json');cases=json.loads(config.read_text())['cases'] if config.exists() else [8,16,32]
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as executor:rows=list(executor.map(inspect,cases))
    atomic('vortex_rhs_report.json',dict(completed=True,rows=rows))


if __name__=='__main__':run()

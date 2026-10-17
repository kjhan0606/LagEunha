"""Compare a numerical RHS with exact Reynolds rates in the same moving cells."""
import json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from adaptive_coupled_mass import atomic,restore
from saved_curved_contact import checked_json
from simple_wave_exact import physical_raw,fields,GAMMA
from simple_wave_evolution import initial,install_domain,score
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


def density_raw(points,rho,velocity,pressure):
    Pi=rho[:,None,None]*velocity[:,:,None]*velocity[:,None,:]
    return np.column_stack([rho,rho[:,None]*points,
        (rho[:,None,None]*points[:,:,None]*points[:,None,:]).reshape(-1,4),rho[:,None]*velocity,
        (rho[:,None,None]*points[:,:,None]*velocity[:,None,:]).reshape(-1,4),Pi.reshape(-1,4),
        pressure/(GAMMA-1)+.5*np.trace(Pi,axis1=1,axis2=2)])


def exact_rates(cells,edges,sites,gen,z,t,amplitude,order=12,field_function=fields,raw_function=physical_raw):
    rates=[];areas=[];area_rates=np.zeros(len(cells))
    for poly in cells:
        raw,points,w,rho,velocity,pressure=raw_function(poly,t,amplitude,order)
        _,_,_,rt,vt,pt=field_function(points,t,amplitude,True)
        dm=w*rt;momentum_time=rt[:,None]*velocity+rho[:,None]*vt
        Pi_time=rt[:,None,None]*velocity[:,:,None]*velocity[:,None,:]+rho[:,None,None]*(
            vt[:,:,None]*velocity[:,None,:]+velocity[:,:,None]*vt[:,None,:])
        energy_time=pt/(GAMMA-1)+.5*rt*np.sum(velocity*velocity,axis=1)+rho*np.sum(velocity*vt,axis=1)
        rates.append(np.r_[dm.sum(),dm@points,(points.T@(dm[:,None]*points)).ravel(),w@momentum_time,
            (points.T@(w[:,None]*momentum_time)).ravel(),np.einsum('n,nij->ij',w,Pi_time).ravel(),w@energy_time])
        areas.append(w.sum())
    rates=np.asarray(rates);g,w=np.polynomial.legendre.leggauss(order);g=(g+1)/2;w=w/2
    for i,j,start,end,n,length in edges:
        points=start+g[:,None]*(end-start)
        if j>=0:
            delta=sites[j]-sites[i]
            wn=(delta@gen[j]-(points-sites[i])@(gen[j]-gen[i])+.5*(z[i]-z[j]))/np.linalg.norm(delta)
        else:wn=np.asarray([model.boundary_velocity(p)@n for p in points])
        rho,v,P=field_function(points,t,amplitude);change=length*((w*wn)@density_raw(points,rho,v,P))
        volume=length*(w@wn);rates[i]+=change;area_rates[i]+=volume
        if j>=0:rates[j]-=change;area_rates[j]-=volume
    return rates,np.asarray(areas),area_rates


def inspect(entry):
    started=time.monotonic();install_domain();n=entry['n'];amp=.02
    if 'initial' in entry:
        source=checked_json(entry['initial'],entry['initial_sha256'])
        mass=np.asarray(source['mass']);state=np.asarray(source['state']);L=[np.asarray(v) for v in source['factors']]
        if entry.get('evaluate_fresh_rhs'):
            dy,diag=kinetic.rhs(0.,state,mass,covariance_factors=L,face_order=4);source_reused=False
        else:
            trial=checked_json(entry['trial'],entry['trial_sha256'])
            record=restore(trial['record']);np.testing.assert_array_equal(record['y'],state)
            dy=record['d0'];source_reused=True
    else:
        mass,state,L,adjustment=initial(n,amp)
        dy,diag=kinetic.rhs(0.,state,mass,covariance_factors=L,face_order=4)
        source_reused=False
    if not source_reused:
        atomic('n'+str(n)+'_initial_state.json',dict(n=n,t=0.,mass=mass,state=state,factors=L,
            physical_initial_condition=dict(kind='nonlinear_simple_wave',amplitude=amp,n=n),source_state_reused='initial' in entry))
        import hashlib
        atomic('n'+str(n)+'_source_rhs.json',dict(t=0.,state_sha256=hashlib.sha256(state.tobytes()).hexdigest(),
            mass_sha256=hashlib.sha256(mass.tobytes()).hexdigest(),dy=dy,diag=diag))
    a=state[:-4].reshape(-1,20);d=dy[:-4].reshape(-1,20)
    cells,edges=model.geometry(a[:,:2],a[:,2],0.)
    reference,areas,area_rates=exact_rates(cells,edges,a[:,:2],d[:,:2],d[:,2],0.,amp,12)
    fine,_,_=exact_rates(cells,edges,a[:,:2],d[:,:2],d[:,2],0.,amp,16)
    numerical=np.column_stack([np.zeros(len(mass)),d[:,3:]])
    P=(GAMMA-1)*(a[:,19]-.5*(a[:,15]+a[:,18]))/areas
    numerical_P=(GAMMA-1)*(d[:,19]-.5*(d[:,15]+d[:,18]))/areas-P*area_rates/areas
    exact_P=(GAMMA-1)*(reference[:,17]-.5*(reference[:,13]+reference[:,16]))/areas-P*area_rates/areas
    exact_velocity=(reference[:,7:9]-a[:,9:11]/mass[:,None]*reference[:,0,None])/mass[:,None]
    velocity=d[:,9:11]/mass[:,None]
    errors=dict(pressure_mean_rate_L1=float(areas@abs(numerical_P-exact_P)/amp),
        velocity_mean_rate_L1=float(areas@np.linalg.norm(velocity-exact_velocity,axis=1)/amp),
        mass_rate_L1=float(abs(reference[:,0]).sum()/amp),
        momentum_rate_L1=float(np.linalg.norm(numerical[:,7:9]-reference[:,7:9],axis=1).sum()/amp),
        energy_rate_L1=float(abs(numerical[:,17]-reference[:,17]).sum()/amp))
    quadrature=float(np.max(np.sum(abs(reference-fine),axis=0))/amp)
    result=dict(n=n,errors_over_amplitude=errors,reference_quadrature_difference=quadrature,
        reference_mass_rate_global=float(abs(reference[:,0].sum())),
        reference_momentum_rate_global=float(np.max(abs(reference[:,7:9].sum(0)))),
        reference_energy_rate_global=float(abs(reference[:,17].sum())),
        source_rhs_reused=source_reused,seconds=time.monotonic()-started,
        note='Exact nonlinear Euler time derivatives plus Reynolds boundary motion, with the actual numerical material generator/weight rates. Reference never enters numerical RHS.')
    if not source_reused:
        result['numerical_closure']={k:v for k,v in diag.items() if k.startswith(('velocity_','pressure_')) or k in ('mass_flux','mass_root_backend')}
    if entry.get('measure_initial_fields'):result['initial_fields']=score(0.,state,mass,amp,L)
    atomic('n'+str(n)+'_rhs_consistency.json',result);print(result,flush=True);return result


def run():
    cfg=json.loads(Path('rhs_consistency_config.json').read_text())
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as executor:
        rows=list(executor.map(inspect,cfg['cases']))
    atomic('rhs_consistency_report.json',dict(completed=True,rows=rows))


if __name__=='__main__':run()

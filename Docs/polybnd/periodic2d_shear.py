"""Steady sinusoidal shear diagnostic using unchanged 2D fixed-mass ALE RHS.
Conservative-average IC and explicitly different resolved-pressure control IC.
No energy adjustment after initialization, no forcing or steering.
"""
import json
from pathlib import Path
import numpy as np
from periodic2d_hllc import initial as base_initial,mesh,quadrature,rhs,GAMMA

P0=1/GAMMA


def shear_integrals(polygons,amplitude):
    means=[]; squares=[]; modes=[]
    for poly in polygons:
        q,w=quadrature(poly); v=amplitude*np.sin(2*np.pi*q[:,1])
        means.append(w@v/w.sum()); squares.append(w@(v*v)/w.sum())
        modes.append(w@np.exp(-2j*np.pi*q[:,1]))
    return np.array(means),np.array(squares),np.array(modes)


def initial(n,amplitude,closure):
    mass,y=base_initial(n,'uniform',np.zeros(2)); polygons,area,_=mesh(y[:,:2],y[:,2])
    mean,square,_=shear_integrals(polygons,amplitude)
    y[:,3]=mass*mean; y[:,4]=0.
    variance=square-mean*mean
    if np.min(variance)<-1e-13: raise RuntimeError('negative quadrature variance')
    if closure=='conservative': kinetic_density=.5*square
    elif closure=='resolved': kinetic_density=.5*mean*mean
    else: raise ValueError(closure)
    y[:,5]=area*(P0/(GAMMA-1)+kinetic_density)
    return mass,y,dict(unresolved_kinetic=float(.5*(area@variance)),
        pressure_bias_max=float(.5*(GAMMA-1)*np.max(variance)),
        pressure_bias_mean=float(.5*(GAMMA-1)*(area@variance)))


def diagnose(y,mass,amplitude):
    polygons,area,edges=mesh(y[:,:2],y[:,2],audit=True)
    mean,square,modes=shear_integrals(polygons,amplitude)
    velocity=y[:,3:5]/mass[:,None]; kinetic=.5*np.sum(y[:,3:5]**2,axis=1)/mass
    pressure=(GAMMA-1)*(y[:,5]-kinetic)/area; pmean=float(area@pressure)
    projected_pressure=P0+.5*(GAMMA-1)*(square-mean*mean)
    rms=lambda a:float(np.sqrt(area@(a*a)))
    mode=np.sum(y[:,3]/area*modes)
    topology=frozenset((i,j) for i,j,*_ in edges)
    return dict(pressure_mean=pmean,pressure_std=rms(pressure-pmean),pressure_rms_from_p0=rms(pressure-P0),
        pressure_rms_from_projected_exact=rms(pressure-projected_pressure),
        projected_pressure_mean=float(area@projected_pressure),projected_pressure_std=rms(projected_pressure-area@projected_pressure),
        pressure_min=float(pressure.min()),pressure_max=float(pressure.max()),density_rms=rms(mass/area-1),
        transverse_velocity_rms=rms(velocity[:,1]),velocity_x_error_rms=rms(velocity[:,0]-mean),
        mode_amplitude=float(2*abs(mode)),mode_phase=float(np.angle(mode)),
        kinetic=float(kinetic.sum()),internal=float(y[:,5].sum()-kinetic.sum()),total_energy=float(y[:,5].sum()),
        momentum=y[:,3:5].sum(0).tolist(),area_min=float(area.min())),topology


def run(n,amplitude,closure,refine=1,end=.15):
    label=f'n{n}_a{amplitude}_{closure}_r{refine}'
    mass,y,bias=initial(n,amplitude,closure); start,topology=diagnose(y,mass,amplitude)
    steps=int(np.ceil(end*n*(1+amplitude)/.12))*refine; dt=end/steps
    stats=dict(mass_flux=0.,pressure_min=float('inf'),area_min=float('inf'),iterations=0)
    telemetry=Path(label+'.history.jsonl'); checkpoint=Path(label+'.checkpoint.json')
    topology_steps=0; edge_additions=0; edge_removals=0
    def observed(state):
        dy,d=rhs(state,mass)
        for k in stats:stats[k]=(min if k in ('pressure_min','area_min') else max)(stats[k],d[k])
        return dy
    def save_checkpoint(step):
        checkpoint.write_text(json.dumps(dict(step=step,time=step*dt,dt=dt,mass=mass.tolist(),state=y.tolist())))
    with telemetry.open('x') as output:
        output.write(json.dumps(dict(step=0,time=0.,added_edges=0,removed_edges=0,**start))+'\n');output.flush()
        save_checkpoint(0)
        for step in range(1,steps+1):
            k1=observed(y); k2=observed(y+dt*k1/2); k3=observed(y+dt*k2/2); k4=observed(y+dt*k3)
            y+=dt*(k1+2*k2+2*k3+k4)/6
            d,current=diagnose(y,mass,amplitude)
            added=len(current-topology); removed=len(topology-current)
            topology_steps+=bool(added or removed);edge_additions+=added;edge_removals+=removed;topology=current
            if abs(d['total_energy']/start['total_energy']-1)>1e-10:raise RuntimeError('energy drift')
            if np.max(abs(np.array(d['momentum'])-start['momentum']))>1e-10:raise RuntimeError('momentum drift')
            output.write(json.dumps(dict(step=step,time=step*dt,added_edges=added,removed_edges=removed,**d))+'\n');output.flush()
            if step%5==0 or step==steps:
                save_checkpoint(step)
                print(f'{label} step {step}/{steps} t={step*dt:.6f} Pstd={d["pressure_std"]:.6g} vyRMS={d["transverse_velocity_rms"]:.6g}',flush=True)
    observed(y)
    result=dict(n=n,amplitude=amplitude,closure=closure,refine=refine,end=end,steps=steps,dt=dt,**stats,
        initial_bias=bias,initial=start,final=d,topology_steps=topology_steps,edge_additions=edge_additions,edge_removals=edge_removals,
        mode_retention=d['mode_amplitude']/start['mode_amplitude'],kinetic_change=d['kinetic']-start['kinetic'],
        internal_change=d['internal']-start['internal'],energy_relative=d['total_energy']/start['total_energy']-1,
        history=str(telemetry),checkpoint=str(checkpoint),final_state=y.tolist())
    assert stats['mass_flux']<3e-13
    return result


if __name__=='__main__':
    cases=[(8,.5,'conservative',1),(8,.5,'resolved',1),(16,.5,'conservative',1),
           (16,.5,'resolved',1),(8,1.5,'conservative',1),(16,.5,'conservative',2)]
    failures=[]
    with open('results.jsonl','x') as output:
        for args in cases:
            try:result=run(*args);result['status']='complete'
            except Exception as error:
                result=dict(args=args,status='failed',error=repr(error));failures.append(result)
            output.write(json.dumps(result)+'\n');output.flush()
    if failures:raise RuntimeError(f'{len(failures)} shear cases failed; see results.jsonl and checkpoints')

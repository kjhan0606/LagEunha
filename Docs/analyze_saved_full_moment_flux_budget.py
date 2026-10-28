"""Signed instantaneous density/pressure budgets from the single saved RHS.

No evolution, parameter sweep, reference forcing or acceptance changes.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from compare_full_moment_checkpoint_errors import load,REGIONS


def read(path):
    with np.load(path) as f:return {k:f[k] for k in f.files}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--capture',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();assert not args.output.exists()
    capture=json.loads((args.capture/'report.json').read_text())
    assert capture['physical_providers']==1 and capture['accepted_time_steps']==0
    root=Path('/gpfs/kjhan/LagForce/Codex_full_disk_continuous_two_dt_to_t1_20261019/coarse/checkpoints/coarse_0005')
    state=read(root/'accepted.npz');mesh=read(root/'geometry/mesh.npz')
    fields,physical,source=load(root/'physics');s=read(args.capture/'stage_compact.npz');flux=read(args.capture/'signed_flux_volume.npz')
    ld=np.longdouble;n=len(state['ids']);ids=state['ids'];mass=state['mass'].astype(ld)
    assert all(np.array_equal(ids,x['ids']) for x in (s,flux,fields,mesh))
    assert s['initial_state_sha256'].item()==hashlib.sha256(np.ascontiguousarray(state['state']).view(np.uint8)).hexdigest()
    a=state['state'][:-4].reshape(n,20).astype(ld);raw=np.column_stack([mass,a[:,3:]])
    st=s['stage'].astype(ld);area=mesh['moments'][:,0].astype(ld);Vdot=flux['Vdot']
    outer=lambda x,y:x[:,:,None]*y[:,None,:]
    q,mom=raw[:,1:3],raw[:,7:9]
    S=raw[:,3:7].reshape(n,2,2)-outer(q,q)/mass[:,None,None]
    C=raw[:,9:13].reshape(n,2,2)-outer(q,mom)/mass[:,None,None]
    B=np.linalg.solve(S.astype(float),C.astype(float)).transpose(0,2,1).astype(ld)
    qdot,pdot=st[:,1:3],st[:,7:9]
    Sdot=st[:,3:7].reshape(n,2,2)-(outer(qdot,q)+outer(q,qdot))/mass[:,None,None]
    Cdot=st[:,9:13].reshape(n,2,2)-(outer(qdot,mom)+outer(q,pdot))/mass[:,None,None]
    u,sigma,vh=np.linalg.svd(state['factors']);scaled=sigma/sigma.max(axis=1)[:,None]
    ratio=scaled[:,:,None]*scaled[:,None,:]/(scaled[:,:,None]**2+scaled[:,None,:]**2)
    K=st[:,22:26].reshape(n,2,2).astype(float)
    loss=(u@((vh@K@vh.transpose(0,2,1))*ratio)@u.transpose(0,2,1)).astype(ld)
    L=state['factors'].astype(ld);R=L@L.transpose(0,2,1)
    Ac=st[:,56:60].reshape(n,2,2)-loss
    Rdot=Ac@R+R@Ac.transpose(0,2,1)+st[:,26:30].reshape(n,2,2)
    BC=B@Cdot
    factorPi=(outer(pdot,mom)+outer(mom,pdot))/mass[:,None,None]+BC+BC.transpose(0,2,1)-B@Sdot@B.transpose(0,2,1)+Rdot
    trace=lambda x:np.trace(x,axis1=1,axis2=2)
    thermal=raw[:,17]-.5*(raw[:,13]+raw[:,16])
    physics=json.loads((Path('/gpfs/kjhan/LagForce/Codex_full_curved_kepler_IC_20261019/report.json')).read_text())['physics']
    scale=(ld(physics['gamma'])-1)/(ld(physics['P0'])*area)
    rates={'cell_mean_density_fractional_rate':-Vdot/area,'relative_pressure_volume':-scale*thermal*Vdot/area,
        'relative_pressure_raw_thermal':scale*(st[:,17]-.5*(st[:,13]+st[:,16])),
        'relative_pressure_factor_minus_raw':-.5*scale*(trace(factorPi)-st[:,13]-st[:,16]),
        'relative_pressure_direct_gravity_thermal':scale*(st[:,47]-.5*(st[:,43]+st[:,46])),
        'advected_internal_energy_rate':flux['advected_internal_energy_rate'],
        'advected_kinetic_energy_rate':flux['advected_kinetic_energy_rate'],
        'pressure_power_rate':flux['pressure_power_rate'],
        'signed_mass_rate':flux['mass_rate'],'native_mass_rate':st[:,0],
        'signed_minus_native_mass':flux['mass_rate']-st[:,0]}
    rates['relative_pressure_factor_total']=rates['relative_pressure_volume']+rates['relative_pressure_raw_thermal']+rates['relative_pressure_factor_minus_raw']
    rates['signed_energy_minus_native_hydro']=flux['advected_internal_energy_rate']+flux['advected_kinetic_energy_rate']+flux['pressure_power_rate']-(st[:,17]-st[:,47])
    def summary(x,mask):
        x=x[mask];w=mass[mask];k=int(np.argmax(abs(x)))
        return dict(RMS=float(np.sqrt(np.sum(w*x*x)/w.sum())),signed_mean=float(np.sum(w*x)/w.sum()),
                    maximum=float(abs(x[k])),maximum_signed=float(x[k]),maximum_id=int(ids[mask][k]))
    regions={name:{k:summary(x,fields[name]) for k,x in rates.items()} for name in REGIONS if fields[name].any()}
    args.output.mkdir();np.savez(args.output/'signed_rates.npz',ids=ids,**rates)
    report=dict(status='completed_single_saved_state_signed_pressure_density_budget',time=.025,
        source=source,capture=capture,regions=regions,physical_providers_added=0,accepted_time_steps=0,
        limits='Instantaneous algebraic tangent, not a new accepted state or finite-step attribution. Density rate is dlog(m/V)/dt for fixed mass, not the derivative of error against a varying actual-cell reference at contacts. Gravity thermal subset is already contained in raw thermal and is not added again. Factor Pi tangent is infinitesimal exponential/QR, not a replacement integrator.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],cavity=regions['cavity'],inner_contact=regions['inner_contact'])),flush=True)


if __name__=='__main__':main()

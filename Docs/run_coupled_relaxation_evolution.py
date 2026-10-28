"""New moving-cell coupled relaxation experiment; not a long Kepler substitute.

Uses the existing factor/midpoint research integrator (not RK4), with its
unvalidated donor residual-transport closure made explicit.
"""
import argparse,hashlib,json,os,time
from pathlib import Path
import numpy as np
from adaptive_moment_batch import setup,diagnostic
from compatible_moment_rhs import rhs
from factor_moment_evolve import rebuild,advance_factors_exponential
from moment_ale_evolve import GAMMA

def energy_parts(y,mass):
    totals=dict(translation=0.,rotation=0.,deformation=0.,residual=0.,thermal=0.,total=0.)
    J=np.array([[0.,-1.],[1.,0.]])
    for m,a in zip(mass,y[:-4].reshape(len(mass),20)):
        q=a[3:5];p=a[9:11];Pi=a[15:19].reshape(2,2)
        S=a[5:9].reshape(2,2)-np.outer(q,q)/m
        C=a[11:15].reshape(2,2)-np.outer(q,p)/m
        G=np.linalg.solve(S,C).T
        ell=C[0,1]-C[1,0];I=np.trace(S);H=G-(ell/I)*J
        R=Pi-np.outer(p,p)/m-G@S@G.T
        for k,v in dict(translation=p@p/(2*m),rotation=ell**2/(2*I),
            deformation=.5*np.trace(H@S@H.T),residual=.5*np.trace(R),
            thermal=a[19]-.5*np.trace(Pi),total=a[19]).items():totals[k]+=float(v)
    totals['decomposition_residual']=sum(v for k,v in totals.items() if k!='total')-totals['total']
    return totals

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=Path('/gpfs/kjhan/LagForce/Codex_coupled_relaxation_evolution_20261018/results'))
parser.add_argument('--dt',type=float,nargs='+',default=[.001,.0005])
parser.add_argument('--end',type=float,default=.01)
parser.add_argument('--pressure',type=float,default=.6)
parser.add_argument('--transport',choices=['donor','kinetic'],default='donor')
parser.add_argument('--quadratic',action='store_true')
args=parser.parse_args()
if not np.isfinite(args.end) or args.end<=0 or any(not np.isfinite(dt) or dt<=0 or abs(round(args.end/dt)*dt-args.end)>1e-12 for dt in args.dt):
    parser.error('Positive finite end/dt and integer step counts required')
if not np.isfinite(args.pressure) or args.pressure<=0:parser.error('Positive finite pressure required')
if args.quadratic and args.transport!='kinetic':parser.error('Quadratic evolution requires kinetic transport')
if args.transport=='kinetic':
    from kinetic_factor_rhs import rhs as kinetic_rhs
    os.environ.update(LAGFORCE_QUADRATIC_VELOCITY=str(int(args.quadratic)),LAGFORCE_QUADRATIC_BUDGET='projected',
        LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area',
        LAGFORCE_QUADRATIC_PRESSURE='0',LAGFORCE_KINETIC_SMOOTHING='0',
        LAGFORCE_VECTOR_KINETIC_ROOT='0',LAGFORCE_PRESSURE_RIEMANN='hllc_davis')
    def rhs(*a,relaxation_alpha=0.,**kw):
        if relaxation_alpha!=0:raise ValueError('No quadratic relaxation closure implemented')
        return kinetic_rhs(*a,**kw)
output=args.output;output.mkdir(exist_ok=False)
results=[]
for alpha in ([0.] if args.transport=='kinetic' else [0.,1.]):
 for dt in args.dt:
    tag=f'alpha{alpha:g}_dt{dt:g}';mass,y=setup('shear',4,0.)
    a=y[:-4].reshape(len(mass),20)
    # Initial density is exactly one; mass equals initial cell area.
    a[:,19]+=(args.pressure-.6)*mass/(GAMMA-1)
    for i,m in enumerate(mass):
        dv=.01*np.array([(-1.)**i,(-1.)**(i//4)])
        q=a[i,3:5];p=a[i,9:11].copy()
        deltaPi=np.outer(dv,p)+np.outer(p,dv)+m*np.outer(dv,dv)
        a[i,9:11]+=m*dv;a[i,11:15]+=np.outer(q,dv).ravel()
        a[i,15:19]+=deltaPi.ravel();a[i,19]+=.5*np.trace(deltaPi)
    factors=[np.zeros((2,2)) for _ in mass];initial=y.copy();y=rebuild(y,mass,factors)
    adjustment=float(abs(y-initial).max());initial=y.copy()
    input_hash=hashlib.sha256(initial.tobytes()+mass.tobytes()).hexdigest()
    start=time.monotonic();heat=0.;minimum_pressure=float('inf');mass_rate=0.;t=0.;calls=0
    history=[dict(time=0.,added_heat=0.,energy=energy_parts(y,mass))]
    result=dict(alpha=alpha,dt=dt,target_time=args.end,initial_pressure=args.pressure,transport=args.transport,quadratic=args.quadratic,initial_hash=input_hash,initial_roundoff_adjustment=adjustment)
    try:
        for step in range(round(args.end/dt)):
            d,diag=rhs(t,y,mass,linear_pressure=True,covariance_factors=factors,return_covariance_terms=True,relaxation_alpha=alpha);calls+=1
            halfL=advance_factors_exponential(factors,diag,dt/2)
            half=rebuild(y+dt*d/2,mass,halfL)
            dh,mid=rhs(t+dt/2,half,mass,linear_pressure=True,covariance_factors=halfL,return_covariance_terms=True,relaxation_alpha=alpha);calls+=1
            nextL=advance_factors_exponential(factors,mid,dt)
            nexty=rebuild(y+dt*dh,mass,nextL)
            minimum_pressure=min(minimum_pressure,diag['min_pressure'],mid['min_pressure'])
            mass_rate=max(mass_rate,diag['mass_flux'],mid['mass_flux'])
            heat+=dt*mid.get('relaxation_heat_rate',0.)
            y,factors=nexty,nextL;t=(step+1)*dt
            if (step+1)%10==0 or step+1==round(args.end/dt):
                history.append(dict(time=t,added_heat=heat,energy=energy_parts(y,mass)))
        final=diagnostic(t,y,mass,initial)
        final['pressure_departure']=max(abs(final['min_pressure']-args.pressure),abs(final['max_pressure']-args.pressure))
        raw=y[:-4].reshape(len(mass),20)
        c=raw[:,3:5]/mass[:,None];v=raw[:,9:11]/mass[:,None]
        ref=np.column_stack((.2+1.1*c[:,1],np.zeros(len(mass))))
        result.update(status='completed',final=final,
          velocity_departure_mass_RMS=float(np.sqrt(np.average(np.sum((v-ref)**2,axis=1),weights=mass))),
          minimum_factor_eigenvalue=min(float(np.linalg.eigvalsh(L@L.T).min()) for L in factors))
        assert final['min_pressure']>0 and final['ledger_error']<1e-9
    except Exception as error:
        result.update(status='failed',error=repr(error))
    result.update(completed_time=t,calls=calls,seconds=time.monotonic()-start,
                  integrated_added_relaxation_heat=heat,min_stage_pressure=minimum_pressure,max_mass_flux=mass_rate)
    np.savez(output/(tag+'.npz'),time=t,state=y,mass=mass,factors=np.asarray(factors),initial=initial)
    (output/(tag+'.energy.json')).write_text(json.dumps(history,indent=2))
    (output/(tag+'.json')).write_text(json.dumps(result,indent=2))
    results.append(result);print(json.dumps(result),flush=True)
(output/'summary.json').write_text(json.dumps(dict(scope=__doc__,cases=results,
    limits=['Small shear patch; no central gravity, historical Kepler interval or production backend validation.',
            'Added heat excludes pressure work and existing hydro/Riemann dissipation.',
            'Positive factors do not prove spatial realizability or Euler closure.']),indent=2))
assert all(r['status']=='completed' for r in results)

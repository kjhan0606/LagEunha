"""Moving continuous-IC moment patch: hydro, quadratic velocity and gravity.

Stage-admissibility experiment, not a replacement for historical long Kepler.
Potential energy is diagnosed from fitted density, not repaired with heat.
"""
import argparse,hashlib,importlib.util,json,os,sys,time,traceback
from pathlib import Path
import numpy as np

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--dt',type=float,required=True)
parser.add_argument('--end',type=float,default=1e-4)
parser.add_argument('--restart',type=Path)
parser.add_argument('--restart-check-only',action='store_true')
parser.add_argument('--active-radius',type=float,choices=(.55,.60,.65),default=.55)
parser.add_argument('--density-closure',choices=('smooth','contact'),default='smooth')
args=parser.parse_args()
contact=args.density_closure=='contact'
if args.restart_check_only and args.restart is None:parser.error('restart checkpoint required')
if not all(np.isfinite(v) and v>0 for v in (args.dt,args.end)) or abs(round(args.end/args.dt)*args.dt-args.end)>1e-14:
    parser.error('positive finite times and integral step count required')
base=Path('/gpfs/kjhan/LagForce')
for path in (base,base/'Codex_coupled_relaxation_evolution_20261018',base/'Codex_quadratic_split_modules_20261018',Path(__file__).resolve().parent):
    sys.path.insert(0,str(path))
import moment_ale_evolve as model
import fixedmass_hllc
from kinetic_factor_rhs import rhs
from factor_moment_evolve import rebuild,advance_factors_sylvester
from quadratic_upwind_face_nodes import QuadraticUpwindNodes
import quadratic_velocity_reconstruction as quad
from moment_gravity_source import source,covariance_source,covariance_generator
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from audit_canonical_pair_history import gravity_coefficient
from run_canonical_av_branches_20261018 import read
if contact:
    import geometry_uncertainty
    from contact_prior_density import ContactPriorDensity,local_contact_priors
    from contact_upwind_face_nodes import ContactUpwindNodes
    from density_quadrature import rule as density_rule

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
geom=load('geom',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
clip=load('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py')
os.environ.update(LAGFORCE_DENSITY_QUADRATURE_ORDER='12',LAGFORCE_QUADRATIC_VELOCITY='1',
    LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',
    LAGFORCE_VELOCITY_STENCIL='area',LAGFORCE_QUADRATIC_PRESSURE='0',LAGFORCE_KINETIC_SMOOTHING='0',
    LAGFORCE_PRESSURE_RIEMANN='exact_euler_star',LAGFORCE_VECTOR_KINETIC_ROOT='0')
model.GAMMA=float(np.float32(1.666667));fixedmass_hllc.GAMMA=model.GAMMA
path=base/'Codex_continuous_initial_rhs_quadratic_split_20261018/rhs_r0.55_f8.npz'
saved=np.load(path);ids=saved['ids'];mass=saved['mass']
frozen=np.linalg.norm(saved['state'][:-4].reshape(-1,20)[:,:2],axis=1)>=args.active_radius
if args.active_radius==.55:assert np.array_equal(frozen,saved['frozen'])
active=~frozen
y=saved['state'].copy();factors=saved['covariance_factors'].copy();initial=y.copy()
start_time=0.;restart_info=None
z=saved['rates'][:-4].reshape(-1,20)[:,2].copy()
if args.restart is not None:
    checkpoint=args.restart.resolve(strict=True)
    parent=json.loads((checkpoint.parent/'manifest.json').read_text())
    if parent.get('density_closure','smooth')!=args.density_closure:raise ValueError('restart density closure mismatch')
    if parent.get('active_radius',.55)!=args.active_radius:raise ValueError('restart active boundary mismatch')
    if parent['sha256']!=hashlib.sha256(path.read_bytes()).hexdigest():raise ValueError('restart IC mismatch')
    with np.load(checkpoint) as restart:
        if not np.array_equal(restart['ids'],ids) or not np.array_equal(restart['mass'],mass):raise ValueError('restart ID/mass mismatch')
        y=restart['state'].copy();factors=restart['factors'].copy();start_time=float(restart['time'])
        if y.shape!=initial.shape or factors.shape!=(len(ids),2,2):raise ValueError('restart shape mismatch')
        if not np.isfinite(y).all() or not np.isfinite(factors).all() or not np.isfinite(start_time):raise ValueError('nonfinite restart')
        if not np.array_equal(y[:-4].reshape(-1,20)[frozen],initial[:-4].reshape(-1,20)[frozen]):raise ValueError('restart reservoir mismatch')
        if 'weight_rate_seed' in restart:
            z=restart['weight_rate_seed'].copy()
            if z.shape!=(len(ids),) or not np.isfinite(z).all():raise ValueError('invalid weight seed')
        seed_stored='weight_rate_seed' in restart
    if start_time<0 or (not args.restart_check_only and args.end<=start_time):raise ValueError('restart must precede target')
    restart_info=dict(checkpoint=str(checkpoint),sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
                      stored_weight_seed=seed_stored,weight_seed_note='Missing legacy seed uses original IC seed only as nonlinear initial guess.')
    if np.min(np.linalg.eigvalsh(factors@factors.transpose(0,2,1)))<=0:raise ValueError('nonpositive restart covariance')
    if not np.allclose(rebuild(y,mass,factors),y,rtol=1e-12,atol=1e-18):raise ValueError('restart Pi/factor mismatch')
remaining_steps=round((args.end-start_time)/args.dt)
if not args.restart_check_only and abs(start_time+remaining_steps*args.dt-args.end)>1e-14:raise ValueError('nonintegral remaining steps')
ic,t0=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t0==0
allx=np.column_stack([ic['x']-3,ic['y']-3])
builder=EmbeddedLaguerreGeometry(allx,ids,clip.clip)
model.geometry=builder;model.boundary_velocity=lambda point:np.zeros(2)
args.output.mkdir(exist_ok=False)
manifest=dict(scope=__doc__,dt=args.dt,end=args.end,cells=len(ids),active_cells=int(active.sum()),
    density_closure=args.density_closure,
    shared_face_velocity=os.environ.get('LAGFORCE_SHARED_FACE_VELOCITY','0')=='1',
    active_radius=args.active_radius,
    start_time=start_time,restart=restart_info,restart_check_only=args.restart_check_only,
    input=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    integrator='explicit midpoint with Sylvester exponential covariance factors; order unverified',
    limits='New continuous IC representation; frozen outer reservoir. Not original point-state restart, RK4, full CUDA/MPI or long-term validation. Potential transport closure remains unresolved.')
if contact:
    names=('contact_prior_density','geometry_uncertainty','density_quadrature','contact_upwind_face_nodes',
           'kinetic_factor_rhs','quadratic_velocity_reconstruction','quadratic_upwind_face_nodes',
           'moment_gravity_source','embedded_laguerre_geometry','factor_moment_evolve','moment_ale_evolve')
    manifest.update(gravity_quadrature='contact parts positive rule12',
        reconstructed_potential_quadrature='contact parts positive rule24; diagnostic only',
        psi_transport_diagnostic=False,
        module_sha256={name:dict(path=str(Path(sys.modules[name].__file__).resolve()),
            sha256=hashlib.sha256(Path(sys.modules[name].__file__).read_bytes()).hexdigest()) for name in names},
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
events=[]
def emit(row):
    events.append(row);print(json.dumps(row),flush=True)
    with (args.output/'events.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')

def decode(y,L):
    a=y[:-4].reshape(-1,20);raw=np.column_stack([mass,a[:,3:]])
    polys,edges=builder(a[:,:2],a[:,2],0)
    if contact:
        priors=local_contact_priors(polys,raw,edges)
        densities=[ContactPriorDensity(p,r[0],r[1:3],r[3:7].reshape(2,2),meta,12) for p,r,meta in zip(polys,raw,priors)]
        states=[model.decode(p,r,density=d) for p,r,d in zip(polys,raw,densities)]
    else:states=[model.decode(p,r) for p,r in zip(polys,raw)]
    for s,l in zip(states,L):s['R']=l@l.T
    quad.prepare(states,polys,edges,mass,L)
    return a,raw,polys,states

last_contact_classification=None
def evaluate(t,y,L,z):
    global last_contact_classification
    start=time.monotonic();a=y[:-4].reshape(-1,20)
    if contact:
        a,raw,polys,states=decode(y,L)
        densities=[s['density'] for s in states];priors=[d.cut_metadata for d in densities]
        classification=np.array([not p['flat'] for p in priors])
        transitions=0 if last_contact_classification is None else int(np.count_nonzero(classification!=last_contact_classification))
        last_contact_classification=classification
        emit(dict(stage_time=t,contact_nonflat=sum(not p['flat'] for p in priors),
            contact_fit_residual=max(d.residual for d in densities),
            contact_max_theta_norm=max(float(np.linalg.norm(d.theta)) for d in densities),
            contact_classification_changes_since_previous_stage=transitions))
    for it in range(10):
        if contact:
            provider=ContactUpwindNodes(a[:,:2],mass,L,z,8,model.boundary_velocity,priors=priors)
            dy,diag=rhs(t,y,mass,covariance_factors=L,frozen_mask=frozen,face_order=8,face_node_provider=provider,density_provider=lambda i,p,r:densities[i])
        else:
            provider=QuadraticUpwindNodes(a[:,:2],mass,L,z,8,model.boundary_velocity)
            dy,diag=rhs(t,y,mass,covariance_factors=L,frozen_mask=frozen,face_order=8,face_node_provider=provider)
        newz=dy[:-4].reshape(-1,20)[:,2].copy();change=float(abs(newz-z).max())
        emit(dict(stage_time=t,iteration=it,weight_rate_change=change,mass_flux=diag['mass_flux']))
        z=newz
        if change<1e-11:break
    else:raise RuntimeError('moving split/weight iteration did not converge')
    if not contact:a,raw,polys,states=decode(y,L)
    rates=dy[:-4].reshape(-1,20)
    gravity_ledger=np.zeros(4)
    for i in np.flatnonzero(active):
        s=states[i]
        if contact:pts=s['density'].points;dm=s['density'].weighted_mass
        else:
            pts,dv=geom.quadrature(polys[i],12);dm=dv*s['density'].value(pts)
        acc=-gravity_coefficient(np.linalg.norm(pts,axis=1))[:,None]*pts
        vel=s['v']+(pts-s['c'])@s['B'].T+quad.correction(s,pts)
        g=source(pts,dm,vel,acc);rates[i,3:]+=g[1:]
        diag['cov_A'][i]=diag['cov_A'][i]+covariance_generator(L[i],covariance_source(raw[i],g))
        gravity_ledger+=np.r_[g[7:9],g[17],g[10]-g[11]]
    dy[-4:]+=gravity_ledger
    assert np.all(rates[frozen]==0)
    emit(dict(stage_time=t,stage_seconds=time.monotonic()-start,min_pressure=diag['min_pressure'],
              gravity_ledger=gravity_ledger.tolist(),generator_velocity_defect=diag['generator_velocity_defect']))
    return dy,diag,z

def diagnostic(t,y,L):
    a,raw,polys,states=decode(y,L)
    potential=0.
    for i in np.flatnonzero(active):
        if contact:
            for part in states[i]['density'].parts:
                pts,dv=density_rule(part,24)
                potential+=float((dv*states[i]['density'].value(pts))@geom.potential(pts))
        else:
            pts,dv=geom.quadrature(polys[i],24)
            potential+=float((dv*states[i]['density'].value(pts))@geom.potential(pts))
    orig=initial[:-4].reshape(-1,20)
    def totals(b):return np.r_[b[:,9:11].sum(axis=0),b[:,19].sum(),(b[:,12]-b[:,13]).sum()]
    return dict(time=t,min_pressure=min(s['P'] for s in states),
        min_factor_eigenvalue=min(float(np.linalg.eigvalsh(l@l.T).min()) for l in L),
        ledger_error=float(abs(totals(a)-totals(orig)-y[-4:]).max()),
        active_hydro_energy=float(a[active,19].sum()),active_reconstructed_potential=potential,
        active_pressure_relative_RMS=float(np.sqrt(np.mean([(states[i]['P']/float(ic['pressure'][0])-1)**2 for i in np.flatnonzero(active)]))),
        max_generator_displacement=float(np.linalg.norm(a[:,:2]-orig[:,:2],axis=1).max()),
        maximum_relative_weight=float(np.ptp(a[:,2])))

def rebuild_diagnostic(before,after,mass,active):
    """Measure the thermal repartition caused solely by Pi reconstruction."""
    a=before[:-4].reshape(-1,20);b=after[:-4].reshape(-1,20)
    delta=b[:,15:19]-a[:,15:19]
    # Ehyd is unchanged by rebuild; this is not an independently added heat source.
    assert np.array_equal(a[:,19],b[:,19])
    thermal_before=a[:,19]-.5*(a[:,15]+a[:,18])
    du=-.5*(delta[:,0]+delta[:,3])
    selected=du[active];weights=mass[active]
    return dict(maximum_Pi_change=float(abs(delta[active]).max()),
        thermal_change_sum=float(selected.sum()),
        thermal_change_absolute_sum=float(abs(selected).sum()),
        specific_thermal_change_mass_RMS=float(np.sqrt(np.average((selected/weights)**2,weights=weights))),
        maximum_relative_thermal_change=(float(np.max(abs(selected/thermal_before[active])))
            if np.all(thermal_before[active]>0) else None),
        minimum_provisional_thermal=float(thermal_before[active].min()),
        scope='Active cells only; Pi reconstruction at unchanged Ehyd; not total physical heating.')

t=start_time;status='failed'
failure_stage='initial_diagnostic';failure_state=y;failure_factors=factors;failure_time=t;failure_seed=z.copy()
np.savez(args.output/'initial.npz',state=y,factors=factors,mass=mass,ids=ids,time=t,weight_rate_seed=z)
try:
    initial_diagnostic=diagnostic(t,y,factors)
    if initial_diagnostic['ledger_error']>1e-9:raise RuntimeError('restart conservation ledger failed')
    emit(dict(initial=initial_diagnostic))
    for step in range(0 if args.restart_check_only else remaining_steps):
        failure_stage='start_rhs';failure_state=y;failure_factors=factors;failure_time=t;failure_seed=z.copy()
        d,diag,z=evaluate(t,y,factors,z)
        halfL=advance_factors_sylvester(factors,diag,args.dt/2)
        half_provisional=y+args.dt*d/2
        half=rebuild(half_provisional,mass,halfL)
        half_rebuild=rebuild_diagnostic(half_provisional,half,mass,active)
        half[:-4].reshape(-1,20)[frozen]=y[:-4].reshape(-1,20)[frozen]
        failure_stage='midpoint_rhs';failure_state=half;failure_factors=halfL;failure_time=t+args.dt/2;failure_seed=z.copy()
        dh,mid,z=evaluate(t+args.dt/2,half,halfL,z)
        nextL=advance_factors_sylvester(factors,mid,args.dt)
        next_provisional=y+args.dt*dh
        nexty=rebuild(next_provisional,mass,nextL)
        next_rebuild=rebuild_diagnostic(next_provisional,nexty,mass,active)
        nexty[:-4].reshape(-1,20)[frozen]=y[:-4].reshape(-1,20)[frozen]
        failure_stage='accepted_candidate_diagnostic';failure_state=nexty;failure_factors=nextL;failure_time=t+args.dt;failure_seed=z.copy()
        result=diagnostic(t+args.dt,nexty,nextL)
        result['rebuild_thermal_partition']=dict(half=half_rebuild,full=next_rebuild)
        if result['ledger_error']>1e-9:raise RuntimeError('moment conservation ledger failed')
        y,factors=nexty,np.asarray(nextL);t=start_time+(step+1)*args.dt
        np.savez(args.output/f'step{step+1:04d}.npz',state=y,factors=factors,mass=mass,ids=ids,time=t,weight_rate_seed=z)
        emit(dict(accepted=result))
    status='restart_validated' if args.restart_check_only else 'completed'
except Exception as exc:
    np.savez(args.output/'failure_stage.npz',state=failure_state,factors=np.asarray(failure_factors),
        mass=mass,ids=ids,time=failure_time,weight_rate_seed=failure_seed,stage=failure_stage)
    emit(dict(status='failed',completed_time=t,failure_stage=failure_stage,failure_stage_time=failure_time,
        failure_checkpoint='failure_stage.npz',error=repr(exc),traceback=traceback.format_exc()))
finally:
    (args.output/'summary.json').write_text(json.dumps(dict(status=status,completed_time=t,events=events),indent=2))
if status not in ('completed','restart_validated'):sys.exit(1)

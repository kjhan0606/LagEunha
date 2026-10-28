"""One actual CUDA-driven moving midpoint step, with common-state shadows.

Python supplies geometry/reconstruction/massroot at each actual CUDA state.
This is not a full native production solver, RK4, MPI or long-term validation.
"""
import argparse,hashlib,importlib.util,json,os,sys,time,traceback
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
modules=base/'Codex_compiled_midpoint_modules_20261018'
sys.path[:0]=[str(modules),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from moment_compiled_backend import Backend
from moment_bridge_audit import stage_reference,audit_stage,factor_reference,audit_factor
from kinetic_face_capture import capture
import moment_ale_evolve as model
import fixedmass_hllc
import quadratic_velocity_reconstruction as quad
import geometry_uncertainty
from kinetic_factor_rhs import rhs,ALPHA
from contact_prior_density import ContactPriorDensity,local_contact_priors
from contact_upwind_face_nodes import ContactUpwindNodes
from moment_gravity_source import source,covariance_source,covariance_generator
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from audit_canonical_pair_history import gravity_coefficient
from run_canonical_av_branches_20261018 import read

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--library',type=Path,default=base/'Codex_compiled_midpoint_build_20261018/liblagmoment.so')
parser.add_argument('--dt',type=float,default=1e-4)
args=parser.parse_args()
if args.dt!=1e-4:parser.error('This approved connection gate has dt=1e-4.')
args.output.mkdir(exist_ok=False)
def emit(row):
    print(json.dumps(row),flush=True)
    with (args.output/'events.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
def sha_array(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
clip=load('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py')
os.environ.update(LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area',LAGFORCE_QUADRATIC_PRESSURE='0',LAGFORCE_KINETIC_SMOOTHING='0',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star',LAGFORCE_VECTOR_KINETIC_ROOT='0',LAGFORCE_SHARED_FACE_VELOCITY='1',LAGFORCE_CAPTURE_PRESSURE_RATES='1')
model.GAMMA=float(np.float32(1.666667));fixedmass_hllc.GAMMA=model.GAMMA
initial_path=base/'Codex_continuous_initial_rhs_quadratic_split_20261018/rhs_r0.55_f8.npz'
saved=np.load(initial_path);mass=saved['mass'];ids=saved['ids'];y0=saved['state'].copy();L0=saved['covariance_factors'].copy()
a0=y0[:-4].reshape(-1,20);n=len(ids);frozen=np.linalg.norm(a0[:,:2],axis=1)>=.60;active=~frozen
z0=saved['rates'][:-4].reshape(-1,20)[:,2].copy()
coarse_initial=np.load(base/'Codex_contact_prior_evolution_dt1e4_20261018/initial.npz')
for key,value in [('state',y0),('factors',L0),('mass',mass),('ids',ids)]:assert np.array_equal(coarse_initial[key],value)
assert np.array_equal(coarse_initial['weight_rate_seed'],z0)
ic,t=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t==0
allx=np.column_stack([ic['x']-3,ic['y']-3]);builder=EmbeddedLaguerreGeometry(allx,ids,clip.clip)
model.geometry=builder;model.boundary_velocity=lambda point:np.zeros(2)
backends={name:Backend(args.library,kind,threads) for name,kind,threads in [('cpu1','cpu',1),('omp4','cpu',4),('cuda','cuda',1)]}
manifest=dict(scope=__doc__,dt=args.dt,cells=n,active_cells=int(active.sum()),active_radius=.60,
    initial_sha256={key:sha_array(value) for key,value in [('state',y0),('factors',L0),('mass',mass),('ids',ids),('weight_seed',z0)]},mask_sha256=sha_array(frozen),
    source=str(initial_path),source_sha256=hashlib.sha256(initial_path.read_bytes()).hexdigest(),library=str(args.library),library_sha256=hashlib.sha256(args.library.read_bytes()).hexdigest(),
    provider='Python geometry/contact density/quadratic velocity/exact-star face split/global anchored massroot/gravity quadrature12; actual CUDA raw/factor states',
    integrator='explicit midpoint + Sylvester exponential factors + Pi rebuild; not production RK4; order unverified',physical_viscosity=0.,Re='infinity',Psi_mode='absent',
    shadow_scope='Same CUDA-selected stage and same CUDA update input, not independent backend trajectories.')
manifest['runner_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
np.savez(args.output/'initial.npz',state=y0,factors=L0,mass=mass,ids=ids,frozen=frozen,weight_rate_seed=z0,time=0.)
progress={}

def decode(y,L):
    a=y[:-4].reshape(-1,20);raw=np.column_stack([mass,a[:,3:]])
    polygons,edges=builder(a[:,:2],a[:,2],0)
    priors=local_contact_priors(polygons,raw,edges)
    densities=[ContactPriorDensity(p,r[0],r[1:3],r[3:7].reshape(2,2),meta,12) for p,r,meta in zip(polygons,raw,priors)]
    states=[model.decode(p,r,density=d) for p,r,d in zip(polygons,raw,densities)]
    for state,Lcell in zip(states,L):state['R']=Lcell@Lcell.T
    quad.prepare(states,polygons,edges,mass,L)
    return a,raw,polygons,edges,priors,densities,states

def evaluate(t,y,L,z,label):
    start=time.monotonic();folder=args.output/label;folder.mkdir(exist_ok=False)
    progress.clear();progress.update(state=y.copy(),factors=L.copy(),seed=z.copy(),label=label,time=t)
    a,raw,polygons,edges,priors,densities,states=decode(y,L)
    emit(dict(stage=label,time=t,input_state_sha256=sha_array(y),factor_sha256=sha_array(L),mask_sha256=sha_array(frozen),nonflat=sum(not p['flat'] for p in priors),theta_max=max(float(np.linalg.norm(d.theta)) for d in densities),fit_residual=max(d.residual for d in densities)))
    for iteration in range(10):
        nodes=ContactUpwindNodes(a[:,:2],mass,L,z,8,model.boundary_velocity,priors=priors)
        dy,diag=rhs(t,y,mass,covariance_factors=L,frozen_mask=frozen,face_order=8,face_node_provider=nodes,density_provider=lambda i,p,r:densities[i],audit_capture_nodes='records')
        newz=dy[:-4].reshape(-1,20)[:,2].copy();change=float(abs(newz-z).max());z=newz
        progress['seed']=z.copy()
        emit(dict(stage=label,time=t,iteration=iteration,weight_rate_change=change,mass_residual=diag['mass_flux']))
        if change<1e-11:break
    else:raise RuntimeError('provider outer did not converge')
    assert diag['mass_flux']<2e-13 and diag['generator_velocity_defect']==0.
    records,values,states,gen,record_z=diag.pop('audit_face_records')
    assert np.array_equal(record_z,z)
    assert np.array_equal(dy[:-4].reshape(-1,20)[:,:2],gen)
    captured=capture(records,values,model.GAMMA,states)
    rows=[[] for _ in range(n)]
    for k,(i,j) in enumerate(captured['owners']):
        rows[i].append(-(k+1))
        if j>=0:rows[j].append(k+1)
    offsets=np.r_[0,np.cumsum([len(row) for row in rows])].astype(np.int64);incidence=np.concatenate(rows).astype(np.int64)
    gn=[];goff=[0];gravity=[];moment_defect=[]
    for i,s in enumerate(states):
        density=s['density'];points=density.points;dm=density.weighted_mass
        velocity=s['v']+(points-s['c'])@s['B'].T+quad.correction(s,points)
        acceleration=-gravity_coefficient(np.linalg.norm(points,axis=1))[:,None]*points
        gn.append(np.column_stack([points,dm,velocity,acceleration]));goff.append(goff[-1]+len(points));gravity.append(source(points,dm,velocity,acceleration))
        samples,_=quad.velocity_samples(s,L[i],mass[i],points,ALPHA)
        pp=dm@velocity;dd=points.T@(dm[:,None]*velocity);pi=np.einsum('n,nki,nkj->ij',dm/4,samples,samples)
        moment_defect.append(np.r_[pp-raw[i,7:9],dd.ravel()-raw[i,9:13],pi.ravel()-raw[i,13:17]])
    gravity=np.array(gravity);moment_defect=np.array(moment_defect)
    assert abs(moment_defect).max()<1e-11
    pv=np.array([s['P']*s['density'].area for s in states])
    data=dict(node_inputs=captured['inputs'],extra=captured['extra'],owners=captured['owners'],cells=np.column_stack([raw,pv]),offsets=offsets,incidence=incidence,frozen=frozen.astype(np.int64),factors=L.copy(),gravity_nodes=np.concatenate(gn),gravity_offsets=np.asarray(goff,dtype=np.int64))
    reference=stage_reference(data,captured,gravity,covariance_source,covariance_generator)
    hydro_original=np.column_stack([diag['mass_rates'],dy[:-4].reshape(-1,20)[:,3:]])
    # Original Python RHS and independent face/volume reconstruction must
    # agree before any backend comparison. No new RHS is evaluated here.
    hydro_difference=float(abs(hydro_original-reference['hydro']).max())
    reference['original_hydro']=hydro_original
    np.savez(folder/'inputs.npz',**data,generator=gen,weight_rates=z,state=y,mass=mass,ids=ids,time=t)
    np.savez(folder/'reference.npz',**reference,moment_defect=moment_defect)
    outputs={name:backend.stage(data) for name,backend in backends.items()}
    for name,out in outputs.items():np.savez(folder/(name+'.npz'),**out)
    comparison=audit_stage(data,reference,outputs)
    report=dict(stage=label,time=t,seconds=time.monotonic()-start,nodes=len(records),gravity_nodes=len(data['gravity_nodes']),maximum_moment_defect=float(abs(moment_defect).max()),hydro_replay_difference=hydro_difference,mass_residual=diag['mass_flux'],outer_change=change,backend_comparison=comparison)
    (folder/'report.json').write_text(json.dumps(report,indent=2));emit(dict(stage_complete=label,time=t,seconds=report['seconds'],mass_residual=diag['mass_flux']))
    return data,outputs,gen,z,report

def update(label,ybase,Lbase,Llinear,data,outputs,gen,z,h):
    folder=args.output/label;folder.mkdir(exist_ok=False)
    a=ybase[:-4].reshape(-1,20);raw=np.column_stack([mass,a[:,3:]]);st=outputs['cuda']['stage']
    provisional=raw+h*st[:,:18];provisional[:,0]=mass;provisional[frozen]=raw[frozen]
    np.savez(folder/'provisional.npz',raw=provisional,Lbase=Lbase,Llinear=Llinear,h=h)
    A=st[:,56:60].reshape(-1,2,2);K=st[:,22:26].reshape(-1,2,2);Q=np.ascontiguousarray(outputs['cuda']['Q'][:,8:].reshape(-1,2));offset=4*data['offsets'];fr=frozen.astype(np.int64)
    inputs=(Lbase,Llinear,A,K,Q,offset,h,provisional,fr)
    reference=factor_reference(*inputs);result={name:b.factor(*inputs) for name,b in backends.items()}
    comparison=audit_factor(inputs,reference,result)
    for name,out in result.items():np.save(folder/(name+'.npy'),out)
    np.save(folder/'python_reference.npy',reference)
    nexty=ybase.copy();aa=nexty[:-4].reshape(-1,20);aa[:,:2]=a[:,:2]+h*gen;aa[:,2]=a[:,2]+h*z;aa[:,3:]=result['cuda'][:,5:22]
    nexty[-4:]=ybase[-4:]+h*outputs['cuda']['ledgers'].reshape(3,4).sum(0)
    aa[frozen]=a[frozen]
    nextL=result['cuda'][:,:4].reshape(-1,2,2).copy();nextL[frozen]=Lbase[frozen]
    assert np.array_equal(result['cuda'][:,4],mass)
    np.savez(folder/'candidate.npz',state=nexty,factors=nextL,mass=mass,ids=ids,frozen=frozen,h=h)
    (folder/'report.json').write_text(json.dumps(comparison,indent=2));emit(dict(update_complete=label,h=h,minimum_thermal=float(result['cuda'][active,51].min())))
    return nexty,nextL

y=y0.copy();L=L0.copy();half=None;halfL=None;candidate=None;candidateL=None
try:
    names=('kinetic_factor_rhs','kinetic_face_capture','moment_compiled_backend','moment_bridge_audit','moment_factor_error_bounds','contact_prior_density','geometry_uncertainty','quadratic_velocity_reconstruction','contact_upwind_face_nodes','moment_gravity_source','embedded_laguerre_geometry')
    manifest['modules']={name:dict(path=str(Path(sys.modules[name].__file__).resolve()),sha256=hashlib.sha256(Path(sys.modules[name].__file__).read_bytes()).hexdigest()) for name in names}
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
    d0,o0,g0,z0,r0=evaluate(0.,y0,L0,z0,'stage0')
    half,halfL=update('half_update',y0,L0,L0,d0,o0,g0,z0,args.dt/2)
    dh,oh,gh,zh,rh=evaluate(args.dt/2,half,halfL,z0,'midpoint')
    candidate,candidateL=update('full_update',y0,L0,halfL,dh,oh,gh,zh,args.dt)
    a,raw,polygons,edges,priors,densities,states=decode(candidate,candidateL)
    assert np.isfinite(candidate).all() and np.isfinite(candidateL).all()
    assert np.array_equal(a[frozen],a0[frozen]) and np.array_equal(candidateL[frozen],L0[frozen])
    assert np.all(raw[:,0]==mass)
    min_R=min(float(np.linalg.eigvalsh(l@l.T).min()) for l in candidateL);assert min_R>0
    thermal=raw[:,17]-.5*(raw[:,13]+raw[:,16]);assert np.all(thermal>0)
    def totals(flat):
        a=flat[:-4].reshape(-1,20);return np.r_[a[:,9:11].sum(0),a[:,19].sum(),(a[:,12]-a[:,13]).sum()]
    ledger_residual=totals(candidate)-totals(y0)-(candidate[-4:]-y0[-4:])
    def absolute_totals(flat):
        a=flat[:-4].reshape(-1,20);return np.r_[abs(a[:,9:11]).sum(0),abs(a[:,19]).sum(),(abs(a[:,12])+abs(a[:,13])).sum()]
    stage_raw=oh['cuda']['stage'][:,:18]
    absolute_rates=np.r_[abs(stage_raw[:,7:9]).sum(0),abs(stage_raw[:,17]).sum(),(abs(stage_raw[:,10])+abs(stage_raw[:,11])).sum()]
    ledger_scale=absolute_totals(candidate)+absolute_totals(y0)+args.dt*(absolute_rates+abs(oh['cuda']['ledgers'].reshape(3,4)).sum(0))
    accepted_bound=(128+4*n)*np.finfo(float).eps*ledger_scale+args.dt*np.asarray(rh['backend_comparison']['cuda']['global_ledger_bound'])
    assert np.all(abs(ledger_residual)<accepted_bound)
    np.savez(args.output/'accepted.npz',state=candidate,factors=candidateL,mass=mass,ids=ids,frozen=frozen,time=args.dt,weight_rate_seed=zh)
    summary=dict(status='completed',accepted_time=args.dt,stage_times=[0.,args.dt/2],minimum_R_eigenvalue=min_R,minimum_thermal=float(thermal.min()),minimum_pressure=min(s['P'] for s in states),ledger_residual=ledger_residual.tolist(),ledger_bound=accepted_bound.tolist(),pressure_relative_RMS=float(np.sqrt(np.mean([(states[i]['P']/float(ic['pressure'][0])-1)**2 for i in np.flatnonzero(active)]))),fit_residual=max(d.residual for d in densities),nonflat=sum(not p['flat'] for p in priors),limits=__doc__)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2));emit(summary)
except Exception as error:
    np.savez(args.output/'failure_state.npz',initial=y0,initial_factors=L0,half=half if half is not None else np.array([]),half_factors=halfL if halfL is not None else np.array([]),candidate=candidate if candidate is not None else np.array([]),candidate_factors=candidateL if candidateL is not None else np.array([]),provider_state=progress.get('state',np.array([])),provider_factors=progress.get('factors',np.array([])),provider_seed=progress.get('seed',np.array([])),mass=mass,ids=ids,frozen=frozen)
    summary=dict(status='failed',error=repr(error),traceback=traceback.format_exc(),provider_stage=progress.get('label'),provider_time=progress.get('time'),limits=__doc__)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2));emit(summary);raise

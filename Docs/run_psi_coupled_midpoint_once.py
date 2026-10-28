"""Actual two-provider Psi-coupled midpoint connection gate, h=1e-4.

No fallback, energy correction, old target replay, or long-time claim.
"""
import argparse,ast,hashlib,importlib.util,json,os,sys,time,traceback
from pathlib import Path
from contextlib import contextmanager
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
modules=base/'Codex_psi_coupled_modules_20261019'
original=base/'Codex_compiled_midpoint_modules_20261018'
sys.path[:0]=[str(modules),str(original),str(base/'Codex_face_state_batch_api_modules_20261018'),str(base/'Codex_moment_flux_batch_modules_20261018'),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
os.environ.update(LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area',LAGFORCE_QUADRATIC_PRESSURE='0',LAGFORCE_KINETIC_SMOOTHING='0',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star',LAGFORCE_VECTOR_KINETIC_ROOT='0',LAGFORCE_SHARED_FACE_VELOCITY='1',LAGFORCE_CAPTURE_PRESSURE_RATES='1',LAGFORCE_BATCH_FACE_STATE='1',LAGFORCE_VECTOR_MOMENT_FLUX='1',LAGFORCE_BATCH_FACE_RECORDS='1',LAGFORCE_NATIVE_EULER_STAR_LIB=str(modules/'libnative_euler_star.so'))
from moment_compiled_backend import Backend
from moment_bridge_audit import stage_reference,audit_stage,factor_reference,audit_factor
from kinetic_face_capture import capture as original_capture
import moment_ale_evolve as model
import fixedmass_hllc
import quadratic_velocity_reconstruction as quad
import geometry_uncertainty
from kinetic_factor_rhs import rhs as original_rhs,ALPHA
from contact_prior_density import ContactPriorDensity,local_contact_priors
from contact_potential_density import ContactPotentialDensity,features
from saved_contact_density_restore import restore_density
from saved_potential_taylor import taylor_coefficients
from contact_upwind_face_nodes import ContactUpwindNodes
from moment_gravity_source import source,covariance_source,covariance_generator
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from audit_canonical_pair_history import gravity_coefficient
from run_canonical_av_branches_20261018 import read
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--library',type=Path,default=base/'Codex_compiled_midpoint_build_20261018/liblagmoment.so');args=parser.parse_args();args.dt=1e-4
args.output.mkdir(exist_ok=False)
assert os.environ['OPENBLAS_NUM_THREADS']==os.environ['OMP_NUM_THREADS']=='1'
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
clip=load('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py')
potential=load('coupled_potential',base/'Codex_weight_remap_moments_centered_20261017/audit.py').potential
model.GAMMA=float(np.float32(1.666667));fixedmass_hllc.GAMMA=model.GAMMA
root=base/'Codex_compiled_midpoint_once_20261018';initial=np.load(root/'initial.npz')
y0=initial['state'].copy();L0=initial['factors'].copy();mass=initial['mass'];ids=initial['ids'];z0=initial['weight_rate_seed'].copy();frozen=initial['frozen'].astype(bool);active=~frozen;n=len(ids);a0=y0[:-4].reshape(-1,20)
ic,t=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t==0
builder=EmbeddedLaguerreGeometry(np.column_stack([ic['x']-3,ic['y']-3]),ids,clip.clip);model.geometry=builder;model.boundary_velocity=lambda point:np.zeros(2)
backends={name:Backend(args.library,kind,threads) for name,kind,threads in [('cpu1','cpu',1),('omp4','cpu',4),('cuda','cuda',1)]}
archive=base/'Codex_numerical_Psi_endpoints_v2_20261018'
initial_fits=json.loads((archive/'initial_fits.json').read_text());initial_priors=json.loads((archive/'initial_priors.json').read_text())
integrals=np.load(base/'Codex_potential_polynomial_20261018/initial_integrals.npz')['integrals']
Psi0=integrals[1,:,0]+integrals[1,:,1];current_Psi=Psi0.copy();current_label='initial';progress={};timers={};provider_calls=0;rhs_rows=[];decoded={}
ld=np.longdouble

def sha_array(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def emit(row):
 print(json.dumps(row),flush=True)
 with (args.output/'events.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
@contextmanager
def region(name):
 begin=time.perf_counter();cpu=time.process_time()
 try:yield
 finally:
  row=timers.setdefault(name,dict(calls=0,wall_seconds=0.,process_cpu_seconds=0.));row['calls']+=1;row['wall_seconds']+=time.perf_counter()-begin;row['process_cpu_seconds']+=time.process_time()-cpu

original_prepare=quad.prepare
original_source=source

def timed_prepare(*a,**kw):
 with region('velocity_prepare'):return original_prepare(*a,**kw)
quad.prepare=timed_prepare

def source(*a,**kw):
 # Gravity source is included in the larger gravity/moment-audit region.
 return original_source(*a,**kw)
def capture(*a,**kw):
 with region('capture'):return original_capture(*a,**kw)
def rhs(*a,**kw):
 with region('rhs_outer'):
  out,diag=original_rhs(*a,**kw)
 rhs_rows.append(dict(label=current_label,face_record=diag['face_record_total_timing'],accumulation=diag['moment_accumulation_timing'],native=diag['face_record_batch'],massroot=diag['coupled_massroot_timing']))
 return out,diag

def decode(y,L):
 a=y[:-4].reshape(-1,20);raw=np.column_stack([mass,a[:,3:]])
 with region('geometry'):polygons,edges=builder(a[:,:2],a[:,2],0)
 with region('prior'):priors=local_contact_priors(polygons,raw,edges)
 folder=args.output/(current_label+'_density');folder.mkdir(exist_ok=False)
 np.savez(folder/'targets.npz',Psi=current_Psi,raw=raw,ids=ids,psi_constrained_mask=active)
 (folder/'priors.json').write_text(json.dumps(priors))
 densities=[];parameters=[]
 for i,(p,r,meta) in enumerate(zip(polygons,raw,priors)):
  with region('baseline_density'):
   if current_label=='initial':
    assert all(meta[k]==initial_priors[i][k] for k in meta if k!='coordinate_error_model')
    old=restore_density(p,r,initial_priors[i],initial_fits[i])
   else:old=ContactPriorDensity(p,r[0],r[1:3],r[3:7].reshape(2,2),meta,12)
  fit=old
  if active[i]:
   with region('potential_density_fit'):
    taylor=taylor_coefficients(old.origin[None,:],potential,gravity_coefficient);s=ld(old.scale);H=taylor['H'][0].astype(ld)
    reference=np.r_[ld(taylor['Phi'][0]),s*taylor['gradient'][0].astype(ld),.5*s*s*H[0,0],s*s*H[0,1],.5*s*s*H[1,1]]
    fit=ContactPotentialDensity(p,r[0],r[1:3],r[3:7].reshape(2,2),meta,old,potential,current_Psi[i],reference,24,(.25,.5))
  row=dict(id=int(ids[i]),active=bool(active[i]),theta=fit.theta.tolist(),logZ=fit.logZ,residual=fit.residual,iterations=fit.iterations,condition=fit.condition,origin=fit.origin.tolist(),scale=fit.scale)
  if active[i]:row.update(reference=[str(v) for v in fit.reference],projection_remainder=[str(v) for v in fit.remainder_projection],sigma=str(fit.sigma),rank_error=str(fit.rank_error),target=str(fit.Psi_target))
  with (folder/'parameters.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
  densities.append(fit);parameters.append(row)
  if (i+1)%500==0:emit(dict(density=current_label,cells=i+1))
 with region('volume_decode'):states=[model.decode(p,r,density=d) for p,r,d in zip(polygons,raw,densities)]
 for state,Lcell in zip(states,L):state['R']=Lcell@Lcell.T
 quad.prepare(states,polygons,edges,mass,L)
 decoded[current_label]=(raw,polygons,densities,states)
 return a,raw,polygons,edges,priors,densities,states

# Preserve tested compiled stage/factor connection maps exactly. Only add one
# large gravity-region timer around the original moment/gravity node loop.
runner=original/'run_compiled_midpoint_once.py';tree=ast.parse(runner.read_text())
functions=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in ('evaluate','update')]
evaluate_fn=next(f for f in functions if f.name=='evaluate')
idx=next(i for i,node in enumerate(evaluate_fn.body) if isinstance(node,ast.For) and isinstance(node.target,ast.Tuple) and ast.unparse(node.target)=='(i, s)')
loop=evaluate_fn.body[idx];evaluate_fn.body[idx]=ast.With(items=[ast.withitem(context_expr=ast.Call(func=ast.Name(id='region',ctx=ast.Load()),args=[ast.Constant('gravity_and_moment_audit')],keywords=[]))],body=[loop])
exec(compile(ast.fix_missing_locations(ast.Module(body=functions,type_ignores=[])),str(runner)+'[preserved maps + region timer]','exec'),globals())

# Independent q48 of all cells; frozen changes are measured, never overwritten.
def independent(label):
 raw,polygons,densities,states=decoded[label];got=[];errors=[]
 from positive_radial_quadrature import rule
 with region('independent_q48'):
  for i,d in enumerate(densities):
   if active[i]:pts,w,prior=d.quadrature(48)
   else:
    parts=[]
    for p in d.parts:
     pts0,w0=rule(p,48,(.25,.5));meta=d.cut_metadata
     rho=meta['value'] if meta['flat'] else (meta['high'] if (np.mean(p,axis=0)@np.asarray(meta['normal'])>meta['offset']) else meta['low'])
     parts.append((pts0,w0,np.full(len(w0),rho)))
    pts=np.concatenate([p[0] for p in parts]);w=np.concatenate([p[1] for p in parts]);prior=np.concatenate([p[2] for p in parts])
   dm=(w*d._value_with_prior(pts,prior)).astype(ld);Phi=potential(pts.astype(ld));A=features(pts,d.origin,d.scale);values=dm@A;psi=dm@Phi
   got.append(psi)
   if active[i]:
    err=values/ld(mass[i])-d.target_A;chi=Phi-A@d.reference-A@d.remainder_projection
    erem=(dm@chi-(d.Psi_target-ld(mass[i])*(d.target_A@d.projection)))/(ld(mass[i])*d.sigma)
    errors.append([i,float(abs(err).max()),float(psi-current_Psi[i]),float(erem)])
    assert np.all(abs(values-ld(mass[i])*d.target_A)<=1e-14+1e-10*abs(ld(mass[i])*d.target_A))
    assert abs(psi-current_Psi[i])<=1e-14+1e-10*abs(current_Psi[i])
 np.savez(args.output/(label+'_q48.npz'),Psi=np.array(got,dtype=ld),errors=np.array(errors),ids=ids,active=active)
 return np.array(got,dtype=ld)

def psi_rate(label,data):
 v=data['node_inputs'];owners=data['owners'];i=owners[:,0];j=owners[:,1];internal=j>=0
 def sum4(a):return ((a[:,0]+a[:,1])+a[:,2])+a[:,3]
 fm=v[:,4]*(sum4(v[:,21:25])+sum4(v[:,25:29]));Phi=potential(v[:,:2]);flux=Phi*fm
 face=-np.bincount(i,weights=flux,minlength=n)+np.bincount(j[internal],weights=flux[internal],minlength=n)
 ref=np.load(args.output/label/'reference.npz');gE=ref['gravity'][:,17];rate=face-gE;rate[frozen]=0.
 velL=v[:,5:13].reshape(-1,4,2);velR=v[:,13:21].reshape(-1,4,2)
 el=.5*np.sum(velL*velL,axis=2)+v[:,29,None];er=.5*np.sum(velR*velR,axis=2)+v[:,30,None]
 eflux=v[:,4]*((v[:,21:25][:,None,:]@el[:,:,None])[:,0,0]+(v[:,25:29][:,None,:]@er[:,:,None])[:,0,0]+v[:,31]*v[:,32])
 right_active=np.zeros(len(j),bool);right_active[internal]=active[j[internal]];direction=active[i].astype(int)-right_active.astype(int)
 bh=float((-direction*eflux).sum());bp=float((-direction*flux).sum())
 combined=float((ref['hydro'][active,17]+gE[active]+rate[active]).sum());res=combined-bh-bp
 ops=8*len(v)+4*n;gamma=ops*np.finfo(float).eps/(1-ops*np.finfo(float).eps);bound=gamma*(abs(eflux).sum()+abs(flux).sum()+abs(gE).sum())
 assert abs(res)<=bound
 np.savez(args.output/label/'Psi_rates.npz',rate=rate,unmasked_rate=face-gE,face=face,gE=gE,node_mass_flux=fm,node_Phi_flux=flux,node_hydro_flux=eflux,boundary_hydro=bh,boundary_Phi=bp,residual=res,bound=bound)
 emit(dict(Psi_ledger=label,residual=res,bound=bound))
 return rate,bh+bp

half=None;halfL=None;candidate=None;candidateL=None
manifest=dict(scope=__doc__,dt=args.dt,Psi_mode='independent transported sidecar; constrained active only',psi_constrained_mask_sha256=sha_array(active),initial_sha256={k:sha_array(initial[k]) for k in initial.files},runner_sha256=sha(__file__),original_runner_sha256=sha(runner),snapshot_sha256=sha(modules/'manifest.json'),Re='infinity',physical_viscosity=0,integrator='explicit midpoint + existing Sylvester exp/QR factors, not RK4',provider_budget=2)
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
np.savez(args.output/'initial.npz',state=y0,factors=L0,mass=mass,ids=ids,frozen=frozen,weight_rate_seed=z0,Psi=Psi0,time=0.,psi_constrained_mask=active)
try:
 current_label='initial';provider_calls+=1;d0,o0,g0,z0,r0=evaluate(0.,y0,L0,z0,'stage0');q0=independent('initial');rate0,_=psi_rate('stage0',d0)
 half,halfL=update('half_update',y0,L0,L0,d0,o0,g0,z0,args.dt/2)
 current_Psi=Psi0+ld(args.dt/2)*rate0.astype(ld);current_Psi[frozen]=Psi0[frozen];current_label='half'
 np.savez(args.output/'half_Psi.npz',Psi=current_Psi,rate0=rate0)
 provider_calls+=1;dh,oh,gh,zh,rh=evaluate(args.dt/2,half,halfL,z0,'midpoint');qh=independent('half');rateh,boundary=psi_rate('midpoint',dh)
 candidate,candidateL=update('full_update',y0,L0,halfL,dh,oh,gh,zh,args.dt)
 current_Psi=Psi0+ld(args.dt)*rateh.astype(ld);current_Psi[frozen]=Psi0[frozen];current_label='final'
 np.savez(args.output/'accepted.npz',state=candidate,factors=candidateL,mass=mass,ids=ids,frozen=frozen,weight_rate_seed=zh,Psi=current_Psi,time=args.dt,psi_constrained_mask=active)
 a,raw,polygons,edges,priors,densities,states=decode(candidate,candidateL);qf=independent('final')
 assert provider_calls==2 and np.isfinite(candidate).all() and np.isfinite(candidateL).all()
 assert np.array_equal(a[frozen].view(np.uint64),a0[frozen].view(np.uint64)) and np.array_equal(candidateL[frozen],L0[frozen]) and np.array_equal(current_Psi[frozen],Psi0[frozen])
 assert np.array_equal(raw[:,0],mass);thermal=raw[:,17]-.5*(raw[:,13]+raw[:,16]);assert np.all(thermal>0)
 assert np.array_equal(raw[:,4],raw[:,5]) and np.array_equal(raw[:,14],raw[:,15])
 minimum_R=min(float(np.linalg.eigvalsh(l@l.T).min()) for l in candidateL);assert minimum_R>0
 hydrodelta=np.sum((raw[active,17]-a0[active,19]).astype(ld));psidelta=np.sum(current_Psi[active]-Psi0[active]);energy=hydrodelta+psidelta-ld(args.dt)*ld(boundary)
 reconstructed=hydrodelta+np.sum(qf[active]-q0[active])-ld(args.dt)*ld(boundary)
 report=dict(status='completed_connection_gate',provider_calls=provider_calls,minimum_R=minimum_R,minimum_thermal=float(thermal.min()),active_transported_energy_residual=float(energy),active_reconstructed_energy_residual=float(reconstructed),active_initial_fit_Psi_error=float(np.sum(q0[active]-Psi0[active])),active_final_fit_Psi_error=float(np.sum(qf[active]-current_Psi[active])),frozen_reconstructed_Psi_change=float(np.sum(qf[frozen]-q0[frozen])),frozen_sidecar_change=float(np.sum(current_Psi[frozen]-Psi0[frozen])),timers=timers,rhs_outer=rhs_rows,physical_comparison='pending separate exact-current-polygon postprocessing',limits=__doc__)
 (args.output/'summary.json').write_text(json.dumps(report,indent=2));emit({k:v for k,v in report.items() if k not in ('rhs_outer','timers')})
except Exception as error:
 np.savez(args.output/'failure_state.npz',initial=y0,initial_factors=L0,half=half if half is not None else [],half_factors=halfL if halfL is not None else [],candidate=candidate if candidate is not None else [],candidate_factors=candidateL if candidateL is not None else [],provider_state=progress.get('state',[]),provider_factors=progress.get('factors',[]),Psi=current_Psi,mass=mass,ids=ids,frozen=frozen)
 report=dict(status='failed',error=repr(error),traceback=traceback.format_exc(),provider_calls=provider_calls,label=current_label,timers=timers,rhs_outer=rhs_rows)
 (args.output/'summary.json').write_text(json.dumps(report,indent=2));emit(dict(status='failed',error=repr(error),label=current_label));raise

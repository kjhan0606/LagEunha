"""One actual moving full-disk stage, reusing its prepared physics.

Budget: one explicitly indexed unique stage of the approved five-stage dt bundle.
No time evolution or automatic fallback. Actual line/arc nodes go to the
existing native exact-star, closed root and existing compiled moment stage.
"""
import argparse,hashlib,importlib.util,json,os,resource,subprocess,sys,time,traceback
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_full_disk_arc_roots_modules_v7_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
os.environ.update(LAGFORCE_NATIVE_EULER_STAR_LIB=str(modules/'libnative_euler_star.so'),LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area')
from curved_disk_cell import DiskCell
from compact_potential_density import CompactPotentialDensity,VolumeCache
sys.path.insert(0,str(base/'Codex_full_disk_fast_modules_v1_20261019'))
from curved_face_nodes_fast import build_geometry as _build_geometry,prepare_states
from same_stage_face_cache import SameStageFaceCache
def build_geometry(*a,**kw):return _build_geometry(*a,stage_cache=static_face_cache,**kw)
from face_record_batch_fast import build as records_build
from closed_root_transport_arrays import solve_from_transport
from closed_weight_root import solve
from explicit_stage_inputs import capture,primitive_energy_momentum
from kinetic_factor_rhs import ALPHA,pressure_gradients
from vector_moment_flux import accumulate
from moment_compiled_backend import Backend
import moment_ale_evolve as model
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--backend',choices=('cuda','cpu'),default='cuda');parser.add_argument('--saved-weight-seed',type=Path);parser.add_argument('--state',type=Path,required=True);parser.add_argument('--geometry',type=Path,required=True);parser.add_argument('--density',type=Path,required=True);parser.add_argument('--velocity-gravity',type=Path,required=True);parser.add_argument('--provider-index',type=int,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False);started=time.perf_counter()
def sha_file(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  while part:=f.read(64*1024**2):h.update(part)
 return h.hexdigest()
def sha_array(a):return hashlib.sha256(np.ascontiguousarray(a).view(np.uint8)).hexdigest()
def emit(row):
 print(json.dumps(row),flush=True)
 with (args.output/'events.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
potential=load('physical_potential',base/'Codex_weight_remap_moments_centered_20261017/audit.py').potential
initial_path=args.state;initial=np.load(initial_path);mass=initial['mass'];ids=initial['ids'];L=initial['factors'];state=initial['state'];Psi=initial['Psi'];N=len(ids);positions=state[:-4].reshape(-1,20)[:,:2];raw=np.column_stack([mass,state[:-4].reshape(-1,20)[:,3:]])
source_root=args.velocity_gravity;source_report=json.loads((source_root/'report.json').read_text());assert source_report['status']=='passed_stage_velocity_gravity_preparation';assert sha_file(initial_path)==source_report['state_sha256']
prepared_file=np.load(source_root/'prepared.npz');prepared={k:prepared_file[k] for k in prepared_file.files}
for key,value in [('state',state),('factors',L),('mass',mass),('Psi',Psi),('ids',ids)]:assert np.array_equal(prepared[key],value)
node_path=source_root/'gravity_nodes.npy';assert sha_file(node_path)==source_report['gravity_input_sha256'];gravity_nodes=np.load(node_path,mmap_mode='r')
fitroot=args.density;fitmanifest=json.loads((fitroot/'manifest.json').read_text());fitreport=json.loads((fitroot/'report.json').read_text());assert fitreport['status']=='passed_full_stage_Psi_density_preparation' and fitreport['state_sha256']==sha_file(initial_path);fitrows=[json.loads(line) for line in (fitroot/'cells.jsonl').read_text().splitlines()];assert [r['id'] for r in fitrows]==ids.tolist()
geom_path=args.geometry/'full_exact_supports.json';native=json.loads(geom_path.read_text());mesh=np.load(args.geometry/'mesh.npz');lookup={int(g):i for i,g in enumerate(ids)}
assert np.array_equal(mesh['ids'],ids) and np.array_equal(mesh['sites'],positions) and np.array_equal(mesh['weights'],state[:-4].reshape(-1,20)[:,2])
assert json.loads((args.geometry/'report.json').read_text())['state_sha256']==sha_file(initial_path)
model.GAMMA=float(np.float32(1.666667));cache=VolumeCache(1);stage_key=sha_file(initial_path)+':'+fitmanifest['source_sha256'][str(geom_path)];cache.begin_stage(stage_key);static_face_cache=SameStageFaceCache(stage_key);cells=[];states=[]
for i,row in enumerate(fitrows):
 cell=DiskCell.from_native(native['cells'][i],2.5);d=CompactPotentialDensity(cell,raw[i],row['prior'],row['parameter'],potential,cache,stage_key,int(ids[i]));s=model.decode(cell,raw[i],density=d);s['R']=L[i]@L[i].T;s['velocity_quadratic']={k:prepared[k][i] for k in ('center','scale','projection','coefficient','standard_coefficient','residual_map')};cells.append(cell);states.append(s)
gen=prepared['generator'];gen_check=np.array([model.velocity(s,p) for s,p in zip(states,positions)]);assert np.array_equal(gen_check,gen)
edges=[(lookup[int(a)],lookup[int(b)]) for a,b in mesh['face_keys']];centers,grad=pressure_gradients(states,edges);priors=[r['prior'] for r in fitrows];frozen=np.zeros(N,bool)
library=base/'Codex_compiled_midpoint_build_20261018/liblagmoment.so';backend=Backend(library,args.backend,1 if args.backend=='cuda' else 4)
manifest=dict(scope=__doc__,provider_calls=1,unique_provider_budget=5,time=float(initial['time']),dt_branches_pending=[.001,.0005],backend=args.backend,initial_sha256=sha_file(initial_path),prepared_sha256=sha_file(source_root/'prepared.npz'),gravity_input_sha256=source_report['gravity_input_sha256'],compiled_library_sha256=sha_file(library),runner_sha256=sha_file(__file__),face_snapshot=str(modules),same_state_reuse='Initial density fits and velocity/gravity preparation exactly reused; no initial refit',Re='infinity',physical_viscosity=0.,implementation_role='same_physical_stage_optimization_replay',new_unique_physical_provider_calls=0,optimization_snapshot=str(base/'Codex_full_disk_fast_modules_v1_20261019'))
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2));np.savez(args.output/'input.npz',state=state,factors=L,mass=mass,ids=ids,Psi=Psi,frozen=frozen)
assert 2<=args.provider_index<=5
manifest['unique_physical_provider_index']=args.provider_index;manifest['attempt']=1
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
outer=[];z=initial['weight_rate_seed'].copy();assert z.shape==(N,) and np.isfinite(z).all()
if args.saved_weight_seed is not None:
 saved=np.load(args.saved_weight_seed)
 for key,value in [('state',state),('factors',L),('mass',mass),('Psi',Psi),('ids',ids)]:assert np.array_equal(saved[key],value)
 z=saved['last_weight_rate_seed'].copy();assert z.shape==(N,) and np.isfinite(z).all()
 manifest['same_state_weight_seed']=dict(path=str(args.saved_weight_seed),sha256=sha_file(args.saved_weight_seed),note='Initial guess only; closed root re-solved on newly split nodes, not accepted without verification')
 manifest['attempt']=2;manifest['unique_physical_provider_index']=1
 (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
last_geometry=None;last_cache=None;last_audit=None;last_values=None;captured=None;outputs=None
try:
 emit(dict(status='provider_started',cells=N,prepared_seconds=time.perf_counter()-started,provider_calls=1))
 for iteration in range(10):
  begin=time.perf_counter();geometry=build_geometry(cells,ids,positions,mass,L,states,gen,z,priors);geometry_seconds=time.perf_counter()-begin;last_geometry=geometry
  (args.output/f'arc_roots_outer{iteration}.json').write_text(json.dumps(geometry['arc_root_audit'],indent=2))
  count=len(geometry['points']);ni=2*int((geometry['owners'][:,1]>=0).sum())+int((geometry['owners'][:,1]<0).sum())
  estimated_gpu_bytes=int(gravity_nodes.nbytes+count*(33+20+18)*8+ni*16*8+N*300*8+count*32)
  estimated_host_bytes=int(gravity_nodes.nbytes+count*2500+ni*16*8+4*1024**3)
  if iteration==0:
   mem={}
   for line in Path('/proc/meminfo').read_text().splitlines():
    k,*v=line.split()
    if k=='MemAvailable:':mem['available_bytes']=int(v[0])*1024
   preflight=dict(nodes=count,incidence=ni,gravity_nodes=len(gravity_nodes),estimated_gpu_bytes=estimated_gpu_bytes,estimated_host_bytes=estimated_host_bytes,host=mem,output_policy='Node hash/count plus compact stage/incoming factor data. No full face node dump; bounded failure witness only.')
   assert estimated_host_bytes<mem['available_bytes']/2
   if args.backend=='cuda':
    result=subprocess.run(['nvidia-smi','--query-gpu=memory.free','--format=csv,noheader,nounits'],check=True,text=True,capture_output=True);free=int(result.stdout.strip().splitlines()[0])*1024**2;preflight['GPU_free_bytes']=free;assert estimated_gpu_bytes<free*.8
   (args.output/'resource_preflight.json').write_text(json.dumps(preflight,indent=2));emit(dict(resource_preflight=preflight))
  begin=time.perf_counter();face_cache=prepare_states(geometry,states,L,mass);last_cache=face_cache;state_seconds=time.perf_counter()-begin
  records,record_stats,audit=records_build(face_cache,[],states,positions,gen,frozen,centers,grad,lambda p:np.zeros(2),model.GAMMA,ALPHA,expose_arrays=True,pack_records=False);last_audit=audit
  begin=time.perf_counter();newz,values,norm,root_diag=solve_from_transport(geometry,audit['transport_arrays'],N);root_seconds=time.perf_counter()-begin;last_values=values;change=float(abs(newz-z).max());z=newz
  row=dict(iteration=iteration,nodes=count,geometry_seconds=geometry_seconds,face_state_seconds=state_seconds,record=record_stats,massroot_seconds=root_seconds,mass_residual=norm,weight_rate_change=change,closed_root=root_diag);row['same_stage_face_cache']=static_face_cache.stats();outer.append(row)
  emit({k:v for k,v in row.items() if k not in ('record','closed_root')})
  if change<1e-11:break
  # No arrays from an obsolete outer are retained by the next iteration.
  del records,face_cache,audit,values,last_cache,last_audit,last_values,geometry,last_geometry
  last_cache=last_audit=last_values=last_geometry=None
 else:raise RuntimeError('actual initial provider outer did not converge')
 assert norm<2e-13
 captured=capture(geometry,audit['transport_arrays'],values,states,model.GAMMA)
 # The independent wall order audit uses the same numerical state and flux
 # API; it performs no additional RHS, geometry update or massroot solve.
 wall_results={}
 for order in (8,16):
  wall_geom=build_geometry(cells,ids,positions,mass,L,states,gen,z,priors,arc_order=order,wall_only=True);wall_cache=prepare_states(wall_geom,states,L,mass)
  wall_records,wall_stats,wall_audit=records_build(wall_cache,[],states,positions,gen,frozen,centers,grad,lambda p:np.zeros(2),model.GAMMA,ALPHA,expose_arrays=True)
  from vector_kinetic_root import record_arrays,evaluate as evaluate_mass
  datawall=record_arrays(wall_records);residual,_,wall_values=evaluate_mass(datawall,N,z,False);wl,wr,wn=wall_values
  # Reuse accumulated-node oracle on this boundary only.
  vector_values=[(a,b,float(c)) for a,b,c in zip(wl,wr,wn)]
  wall_dr,wall_ledger,F,H,K,incoming,area=accumulate(wall_records,vector_values,states,mass,model.GAMMA)
  covariance=np.array([sum((a@a.T for a in blocks),np.zeros((2,2))) for blocks in incoming])
  wall_results[order]=dict(dr=wall_dr,ledger=wall_ledger,F=F,H=H,K=K,incoming=covariance,area=area)
  if order==8:
   wall_capture=capture(wall_geom,wall_audit['transport_arrays'],wall_values,states,model.GAMMA);wf,we,wp,wa=primitive_energy_momentum(wall_capture['node_inputs']);Phi=np.asarray(potential(wall_geom['points'].astype(np.longdouble)),dtype=float)
   transport=wall_audit['transport_arrays'];vl=transport['vl'];vr=transport['vr'];ds=wall_geom['record_scale'];p=wall_geom['points'];rnorm=np.linalg.norm(p,axis=1)
   # Cartesian primitive envelopes cover cancellation in cold tangent dot
   # products, mirrored energies, and the native exact-star root tolerance.
   speedL=np.sum(abs(vl),axis=2);speedR=np.sum(abs(vr),axis=2);eL=abs(wall_capture['node_inputs'][:,29]);eR=abs(wall_capture['node_inputs'][:,30])
   cartM=transport['rl']*speedL.sum(1)+transport['rr']*speedR.sum(1);eps=np.finfo(float).eps
   mass_bound=512*eps*ds*cartM
   energy_primitive=transport['rl']*np.sum(speedL*(speedL**2+eL[:,None]),axis=1)+transport['rr']*np.sum(speedR*(speedR**2+eR[:,None]),axis=1)+abs(transport['ps'])*(np.sum(abs(wall_cache['means'][:,0]),axis=1)+np.sqrt(model.GAMMA*transport['pl']/transport['rl']))
   energy_bound=512*eps*ds*energy_primitive+ds*abs(transport['ps'])*wall_audit['speed_error_bound']
   torque=p[:,0]*wp[:,1]-p[:,1]*wp[:,0];torque_bound=512*eps*ds*rnorm*(transport['rl']*np.sum(speedL**2,axis=1)+transport['rr']*np.sum(speedR**2,axis=1)+abs(transport['ps']))
   assert np.all(abs(wf)<=mass_bound) and np.all(abs(we)<=energy_bound) and np.all(abs(torque)<=torque_bound)
   wall_zero=dict(nodes=len(wf),mass_flux_max=float(abs(wf).max(initial=0.)),energy_flux_max=float(abs(we).max(initial=0.)),Phi_flux_max=float(abs(Phi*wf).max(initial=0.)),origin_torque_max=float(abs(torque).max(initial=0.)),mass_bound_sum=float(mass_bound.sum()),energy_bound_sum=float(energy_bound.sum()),torque_bound_sum=float(torque_bound.sum()))
 wall_order={}
 for name in ('dr','F','H','K','incoming','area'):
  a=wall_results[8][name];b=wall_results[16][name];error=abs(a-b);bound=1e-14+1e-10*abs(b);assert np.all(error<=bound),(name,float(error.max()))
  wall_order[name]=dict(maximum_difference=float(error.max(initial=0.)),maximum_bound_fraction=float(np.max(error/bound,initial=0.)))
 data=dict(captured,cells=np.column_stack([raw,np.array([s['P']*s['density'].area for s in states])]),frozen=np.zeros(N,dtype=np.int64),factors=L,gravity_nodes=gravity_nodes,gravity_offsets=prepared['gravity_offsets'])
 begin=time.perf_counter();outputs=backend.stage(data);compiled_seconds=time.perf_counter()-begin;stage=outputs['stage'];assert np.isfinite(stage).all() and not stage[:,93].any()
 gc=np.diff(data['gravity_offsets']);gravity_bound=(64+2*gc[:,None])*np.finfo(float).eps*prepared['gravity_absolute_sum'];gravity_error=abs(stage[:,30:48]-prepared['gravity']);assert np.all(gravity_error<=gravity_bound)
 assert float(abs(stage[:,0]).max())<2e-13
 fm,energy_flux,fp,absolute_ledger=primitive_energy_momentum(captured['node_inputs']);owners=captured['owners'];i=owners[:,0];j=owners[:,1];inside=j>=0
 Phi=np.asarray(potential(captured['node_inputs'][:,:2].astype(np.longdouble)),dtype=float);flux=Phi*fm
 facePsi=-np.bincount(i,weights=flux,minlength=N)+np.bincount(j[inside],weights=flux[inside],minlength=N);Psi_rate=facePsi-stage[:,47]
 ledgers=outputs['ledgers'].reshape(3,4);global_rate=np.r_[stage[:,7:9].sum(0),stage[:,17].sum(),(stage[:,10]-stage[:,11]).sum()];ledger_residual=global_rate-ledgers.sum(0)
 eps=np.finfo(float).eps;ops=128+8*len(fm)+4*N;gamma=ops*eps/(1-ops*eps)
 scale=absolute_ledger+np.r_[abs(stage[:,37:39]).sum(0),abs(stage[:,47]).sum(),(abs(stage[:,40])+abs(stage[:,41])).sum()];scale[3]+=2*abs(raw[:,[14,15]]).sum()
 ledger_bound=gamma*scale;assert np.all(abs(ledger_residual)<=ledger_bound)
 combined=float((stage[:,17]+Psi_rate).sum());boundary=float(ledgers[0,2]-flux[~inside].sum());energy_residual=combined-boundary;energy_bound=gamma*(absolute_ledger[2]+abs(flux).sum()+abs(stage[:,47]).sum());assert abs(energy_residual)<=energy_bound
 gravity_torque=float((stage[:,40]-stage[:,41]).sum());gravity_torque_bound=float(gravity_bound[:,[10,11]].sum()+(128+4*N)*eps*(abs(stage[:,40])+abs(stage[:,41])).sum());assert abs(gravity_torque)<=gravity_torque_bound
 np.savez(args.output/'stage0.npz',stage=stage,incoming_Q=np.ascontiguousarray(outputs['Q'][:,8:].reshape(-1,2)),incoming_offsets=4*captured['offsets'],generator=gen,weight_rates=z,Psi_rates=Psi_rate,ledgers=outputs['ledgers'],mass=mass,ids=ids,initial_state_sha256=sha_array(state))
 report=dict(status='passed_actual_full_moving_stage',scope=__doc__,provider_calls=1,unique_provider_budget=5,time=float(initial['time']),cells=N,nodes=len(fm),gravity_nodes=len(gravity_nodes),outer=outer,final_mass_residual=norm,compiled_mass_rate_max=float(abs(stage[:,0]).max()),generator_velocity_defect=0.,wall_zero=wall_zero,wall_order8_vs16=wall_order,central_gravity_origin_torque=gravity_torque,central_gravity_torque_bound=gravity_torque_bound,global_ledger_residual=ledger_residual.tolist(),global_ledger_bound=ledger_bound.tolist(),combined_hydro_Psi_boundary_rate_residual=energy_residual,combined_energy_bound=energy_bound,same_stage_face_cache=static_face_cache.stats(),compiled_seconds=compiled_seconds,total_seconds=time.perf_counter()-started,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,node_inputs_sha256=sha_array(captured['node_inputs']),node_extra_sha256=sha_array(captured['extra']),gravity_maximum_difference=float(gravity_error.max()),gravity_maximum_bound_fraction=float(np.divide(gravity_error,gravity_bound,out=np.zeros_like(gravity_error),where=gravity_bound>0).max()),limits='One moving stage only. No accepted time step, two-dt comparison or long-time stability claim. Existing compiled moment/factor ABI retained; all fluid active, no physical viscosity, fallback or energy correction.')
except Exception as exc:
 np.savez(args.output/'failure_state.npz',state=state,factors=L,Psi=Psi,mass=mass,ids=ids,last_weight_rate_seed=z)
 if captured is not None:np.savez(args.output/'bounded_node_witness.npz',node_inputs=captured['node_inputs'][:64],extra=captured['extra'][:64],owners=captured['owners'][:64])
 report=dict(status='failed_actual_full_moving_stage',provider_calls=1,error=repr(exc),traceback=traceback.format_exc(),outer=outer,seconds=time.perf_counter()-started)
 (args.output/'report.json').write_text(json.dumps(report,indent=2));emit(dict(status=report['status'],error=repr(exc)));raise
(args.output/'report.json').write_text(json.dumps(report,indent=2));emit({k:v for k,v in report.items() if k not in ('outer','wall_order8_vs16','scope')})

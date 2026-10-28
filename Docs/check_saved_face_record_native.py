"""Actual saved stage0/half Riemann states; no fit, RHS/root or evolution."""
import argparse,hashlib,importlib.util,json,os,sys,time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_face_record_native_modules_v2_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_face_state_batch_api_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
os.environ.update(LAGFORCE_NATIVE_EULER_STAR_LIB=str(modules/'libnative_euler_star.so'),LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area')
import moment_ale_evolve as model
from saved_contact_density_restore import restore_density
from quadratic_velocity_reconstruction import prepare,velocity_samples_and_mean_batch
from kinetic_factor_rhs import ALPHA,pressure_gradients
from face_record_batch import build
from euler_star_pressure import star
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from run_canonical_av_branches_20261018 import read
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
model.GAMMA=float(np.float32(1.666667));model.boundary_velocity=lambda p:np.zeros(2)
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
clip=load('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py')
ic,t0=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t0==0.
reports={};passed=True
manifest=dict(scope=__doc__,newton_fits=0,provider_calls=0,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),modules/'face_record_batch.py',modules/'native_euler_star.py',modules/'native_euler_star.cpp',modules/'libnative_euler_star.so')})
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
for label in ('stage0','midpoint'):
 started=time.perf_counter();folder=base/'Codex_adjacency_prior_impact_20261018'/label;common=np.load(folder/'common_inputs.npz')
 saved=np.load(base/'Codex_compiled_midpoint_once_20261018'/label/'inputs.npz');priors=json.loads((folder/'priors.json').read_text())['old'];fits=json.loads((folder/'fit_residuals.json').read_text())
 raw=common['raw'];a=saved['state'][:-4].reshape(-1,20);mass=saved['mass'];L=saved['factors'];ids=saved['ids'];frozen=saved['frozen'].astype(bool)
 builder=EmbeddedLaguerreGeometry(np.column_stack([ic['x']-3,ic['y']-3]),ids,clip.clip);polys,edges=builder(a[:,:2],a[:,2],float(saved['time']))
 off=common['polygon_offsets'];vertices=common['polygon_vertices'];assert all(np.array_equal(p,vertices[lo:hi]) for p,lo,hi in zip(polys,off[:-1],off[1:]))
 densities=[restore_density(p,r,m,f['old']) for p,r,m,f in zip(polys,raw,priors,fits)]
 states=[model.decode(p,r,density=d) for p,r,d in zip(polys,raw,densities)]
 for s,l in zip(states,L):s['R']=l@l.T
 prepare(states,polys,edges,mass,L);centers,grad=pressure_gradients(states,edges)
 nodes=saved['node_inputs'];owners=saved['owners'];point=nodes[:,:2];normal=nodes[:,2:4];count=len(nodes)
 rho=np.full((count,2),np.nan);samples=np.full((count,2,4,2),np.nan);coords=np.full_like(samples,np.nan);means=np.full((count,2,2),np.nan)
 entries=[[] for _ in ids]
 for k,(i,j) in enumerate(owners):
  entries[i].append((k,0))
  if j>=0:entries[j].append((k,1))
 for cell,items in enumerate(entries):
  idx=np.asarray(items);rows=idx[:,0];sides=idx[:,1];points=point[rows];inward=(2*sides-1)[:,None]*normal[rows]
  rho[rows,sides]=densities[cell].value_trace_batch(points,inward)
  v,xi,mean=velocity_samples_and_mean_batch(states[cell],L[cell],mass[cell],points,ALPHA)
  samples[rows,sides]=v;coords[rows,sides]=xi;means[rows,sides]=mean
 distance_by_pair={}
 for i,j in owners:
  key=(int(i),int(j))
  if key not in distance_by_pair:distance_by_pair[key]=np.linalg.norm(a[j,:2]-a[i,:2]) if j>=0 else 0.
 cache=dict(points=point,normals=normal,owners=owners,rho=rho,samples=samples,coordinates=coords,means=means,
  record_scale=nodes[:,4],record_distance=np.array([distance_by_pair[(int(i),int(j))] for i,j in owners]))
 records,stats,audit=build(cache,edges,states,a[:,:2],saved['generator'],frozen,centers,grad,model.boundary_velocity,model.GAMMA,ALPHA)
 inp=audit['inputs'];skip=audit['skip'];native=audit['star'];scalar=np.zeros((count,2))
 for k in np.flatnonzero(~skip):
  v=inp[k];scalar[k]=star(v[0],[v[1],0.],v[2],v[3],[v[4],0.],v[5],[1.,0.],v[6])
 psm=np.column_stack([np.array([r[13] for r in records]),np.array([r[14] for r in records])])
 comparisons={}
 for name,new,old in (('native_vs_scalar_star',psm,scalar),('native_vs_original_capture',psm,nodes[:,31:33]),
   ('left_internal_energy',inp[:,2]/((model.GAMMA-1)*inp[:,0]),nodes[:,29]),('right_internal_energy',inp[:,5]/((model.GAMMA-1)*inp[:,3]),nodes[:,30]),
   ('samples_L',np.array([r[5] for r in records]),nodes[:,5:13].reshape(-1,4,2)),('samples_R',np.array([r[6] for r in records]),nodes[:,13:21].reshape(-1,4,2))):
  delta=abs(new-old);bound=1e-14+1e-10*abs(old);ok=bool(np.all(delta<=bound));passed &= ok
  comparisons[name]=dict(passed=ok,bitwise_equal=new.tobytes()==old.tobytes(),maximum_absolute_difference=float(delta.max(initial=0)),maximum_tolerance_fraction=float((delta/bound).max(initial=0)))
 envelope=native[~skip,7:9]+8*np.finfo(float).eps*abs(scalar[~skip]);delta=abs(psm[~skip]-scalar[~skip]);envelope_ok=bool(np.all(delta<=envelope));passed &= envelope_ok
 comparisons['root_error_envelope']=dict(passed=envelope_ok,maximum_fraction=float(np.divide(delta,envelope,out=np.zeros_like(delta),where=envelope>0).max(initial=0)))
 np.savez(args.output/(label+'.npz'),inputs=inp,native=native,scalar=scalar,skip=skip,owners=owners,old_star=nodes[:,31:33])
 reports[label]=dict(nodes=count,actual_star_calls=int((~skip).sum()),seconds=time.perf_counter()-started,comparisons=comparisons,native_record_stats=stats)
 print(json.dumps(dict(stage=label,report=reports[label])),flush=True)
report=dict(status='passed' if passed else 'numerical_gate_failed',manifest=manifest,stages=reports,
 limits='Stored original points/normals/measures with exact old polynomial fit restoration. No new face splitting or provider call. Actual full provider comparison is still required. Root envelopes assume ordinary libm rounding, not certified intervals.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));assert passed,'Saved failed result; no fallback'

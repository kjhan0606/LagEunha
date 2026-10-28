"""Diagnose failed saved initial wall8/16 gate; no massroot/RHS/evolution."""
import argparse,hashlib,importlib.util,json,os,sys,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_full_disk_face_modules_v2_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
os.environ['LAGFORCE_NATIVE_EULER_STAR_LIB']=str(modules/'libnative_euler_star.so')
from curved_disk_cell import DiskCell
from compact_potential_density import CompactPotentialDensity,VolumeCache
from curved_face_nodes import build_geometry,prepare_states
from face_record_batch import build
from kinetic_factor_rhs import ALPHA,pressure_gradients
from vector_kinetic_root import record_arrays,evaluate
from vector_moment_flux import accumulate
from quadratic_velocity_reconstruction import velocity_samples_and_mean_batch
import moment_ale_evolve as model
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False);started=time.perf_counter()
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
potential=load('physical_potential',base/'Codex_weight_remap_moments_centered_20261017/audit.py').potential
failure=base/'Codex_full_disk_initial_stage0_20261019';saved=np.load(failure/'failure_state.npz');mass=saved['mass'];ids=saved['ids'];L=saved['factors'];N=len(ids);positions=saved['state'][:-4].reshape(-1,20)[:,:2];raw=np.column_stack([mass,saved['state'][:-4].reshape(-1,20)[:,3:]]);z=saved['last_weight_rate_seed']
pz=np.load(base/'Codex_full_disk_initial_velocity_gravity_20261019/prepared.npz');pd={k:pz[k] for k in pz.files};gen=pd['generator'];fitroot=base/'Codex_full_disk_initial_density_20261019';rows=[json.loads(line) for line in (fitroot/'cells.jsonl').read_text().splitlines()];priors=[r['prior'] for r in rows]
mesh=np.load(base/'Codex_full_curved_disk_mesh_20261019/mesh.npz');wall=np.flatnonzero(mesh['arc_counts']);lookup={int(g):i for i,g in enumerate(ids)};native=json.loads((base/'Codex_circular_sparse_support_full_20261019/full_exact_supports.json').read_text())
model.GAMMA=float(np.float32(1.666667));states=[]
for i,row in enumerate(rows):
 d=SimpleNamespace(mass=float(mass[i]),area=row['parameter']['area'],origin=np.array(row['parameter']['origin']));s=model.decode(None,raw[i],density=d);s['R']=L[i]@L[i].T;s['velocity_quadratic']={k:pd[k][i] for k in ('center','scale','projection','coefficient','standard_coefficient','residual_map')};states.append(s)
cache=VolumeCache(1);cache.begin_stage('saved-wall-diagnostic');cells=[SimpleNamespace(arcs=[],lines=[]) for _ in ids]
for i in wall:
 c=DiskCell.from_native(native['cells'][i],2.5);cells[i]=c;states[i]['density']=CompactPotentialDensity(c,raw[i],rows[i]['prior'],rows[i]['parameter'],potential,cache,'saved-wall-diagnostic',int(ids[i]))
centers,grad=pressure_gradients(states,[(lookup[int(a)],lookup[int(b)]) for a,b in mesh['face_keys']]);result={}
for order in (8,16):
 geometry=build_geometry(cells,ids,positions,mass,L,states,gen,z,priors,arc_order=order,wall_only=True);facecache=prepare_states(geometry,states,L,mass)
 records,stats,audit=build(facecache,[],states,positions,gen,np.zeros(N,bool),centers,grad,lambda p:np.zeros(2),model.GAMMA,ALPHA,expose_arrays=True)
 residual,_,values=evaluate(record_arrays(records),N,z,False);wl,wr,wn=values
 dr,ledger,F,H,K,incoming,area=accumulate(records,[(a,b,float(c)) for a,b,c in zip(wl,wr,wn)],states,mass,model.GAMMA)
 covariance=np.array([sum((a@a.T for a in blocks),np.zeros((2,2))) for blocks in incoming]);result[order]=dict(dr=dr,F=F,H=H,K=K,incoming=covariance,area=area)
 np.savez(args.output/f'order{order}.npz',ids=ids,wall=wall,**result[order])
metrics={};witness=set()
for name in result[8]:
 a=result[8][name];b=result[16][name];error=abs(a-b);bound=1e-14+1e-10*abs(b);cellerror=error.reshape(N,-1).max(1);i=int(np.argmax(cellerror));bad=np.flatnonzero(np.any((error>bound).reshape(N,-1),axis=1))
 metrics[name]=dict(maximum_difference=float(cellerror[i]),maximum_id=int(ids[i]),failed_cells=int(len(bad)),maximum_bound_fraction=float((error/bound).max()),value8=np.asarray(a[i]).tolist(),value16=np.asarray(b[i]).tolist());witness.add(i)
 for k in bad[:8]:witness.add(int(k))
samples=[]
for i in sorted(witness):
 for k,arc in enumerate(cells[i].arcs):
  theta=np.linspace(arc.lo,arc.hi,257,dtype=np.longdouble);n=np.column_stack([np.cos(theta),np.sin(theta)]);points=np.asarray(np.longdouble(arc.radius)*n,dtype=float);normals=np.asarray(n,dtype=float)
  v,xi,mean=velocity_samples_and_mean_batch(states[i],L[i],mass[i],points,ALPHA);speed=(v@normals[:,:,None])[:,:,0];changes=np.sum(speed[:-1]*speed[1:]<0,axis=0)
  row=dict(id=int(ids[i]),arc=k,lo=arc.lo,hi=arc.hi,normal_speed_min=speed.min(0).tolist(),normal_speed_max=speed.max(0).tolist(),sampled_sign_changes=changes.tolist(),K8=result[8]['K'][i].tolist(),K16=result[16]['K'][i].tolist());samples.append(row)
  np.savez(args.output/f'speed_id{int(ids[i])}_arc{k}.npz',theta=theta,normal_speed=speed,coordinates=xi,points=points,normals=normals)
report=dict(status='completed_saved_wall_diagnostic',scope=__doc__,provider_calls=0,evolution_steps=0,massroot_solves=0,metrics=metrics,witness=samples,seconds=time.perf_counter()-started,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),failure/'report.json',failure/'failure_state.npz',modules/'curved_face_nodes.py']},limits='Reevaluates only saved initial wall quadrature at the already specified8/16 orders, plus direct normal-speed sampling. No automatic order increase, new state, closed-root solve or fallback.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(status=report['status'],metrics=metrics,witness=samples),indent=2))

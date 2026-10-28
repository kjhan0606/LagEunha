"""Actual initial density/velocity/positive gravity-node preparation for stage0.

Uses the saved initial Psi fits without refitting. One necessary 2.32 GB
compiled-kernel input is sized before creation, reused by actual stage0.
No face/root/RHS/evolution is evaluated by this preparation alone.
"""
import argparse,hashlib,importlib.util,json,multiprocessing as mp,os,resource,sys,time,traceback
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.spatial import cKDTree
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_full_disk_density_modules_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
os.environ.update(LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area')
from curved_disk_cell import DiskCell
from compact_potential_density import CompactPotentialDensity,VolumeCache
from quadratic_cell_prepare import prepare_cell
import quadratic_velocity_reconstruction as quad
import moment_ale_evolve as model
from kinetic_factor_rhs import ALPHA
from moment_gravity_source import source
from audit_canonical_pair_history import gravity_coefficient
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False);started=time.perf_counter()
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
potential=load('physical_potential',base/'Codex_weight_remap_moments_centered_20261017/audit.py').potential
initial_path=base/'Codex_full_curved_kepler_IC_20261019/initial.npz';initial=np.load(initial_path);mass=initial['mass'];ids=initial['ids'];factors=initial['factors'];raw=np.column_stack([mass,initial['state'][:-4].reshape(-1,20)[:,3:]]);positions=initial['state'][:-4].reshape(-1,20)[:,:2];N=len(ids)
fitroot=base/'Codex_full_disk_initial_density_20261019';fitreport=json.loads((fitroot/'report.json').read_text());fitmanifest=json.loads((fitroot/'manifest.json').read_text());assert fitreport['status']=='passed_full_initial_Psi_density_preparation'
assert hashlib.sha256(initial_path.read_bytes()).hexdigest()==fitreport['initial_sha256']
fit_rows=[json.loads(line) for line in (fitroot/'cells.jsonl').read_text().splitlines()];assert [r['id'] for r in fit_rows]==ids.tolist() and all(r['status']=='passed' for r in fit_rows)
native_path=base/'Codex_circular_sparse_support_full_20261019/full_exact_supports.json';native=json.loads(native_path.read_text());cells=native['cells'];R=2.5
model.GAMMA=float(np.float32(1.666667));states=[model.decode(None,r,density=SimpleNamespace(mass=float(r[0]),area=float(row['parameter']['area']),origin=np.array(row['parameter']['origin']))) for r,row in zip(raw,fit_rows)]
for s,L in zip(states,factors):s['R']=L@L.T
centers=np.array([s['c'] for s in states]);means=np.array([s['v'] for s in states]);second=raw[:,3:7].reshape(N,2,2)/mass[:,None,None];tree=cKDTree(centers)
offsets=np.r_[0,np.cumsum([r['nodes24'] for r in fit_rows])].astype(np.int64);ng=int(offsets[-1]);expected_bytes=ng*7*8
assert os.environ['OPENBLAS_NUM_THREADS']==os.environ['OMP_NUM_THREADS']=='1' and os.environ['OMP_PROC_BIND']=='false'
st=os.statvfs(base);memory={}
for line in Path('/proc/meminfo').read_text().splitlines():
 k,*v=line.split()
 if k in ('MemTotal:','MemAvailable:'):memory[k[:-1]]=int(v[0])*1024
preflight=dict(cells=N,gravity_nodes=ng,required_kernel_input_bytes=expected_bytes,other_output_estimate_bytes=N*2000,free_filesystem_bytes=st.f_bavail*st.f_frsize,memory=memory,workers=4,worker_BLAS_OMP_threads=1,allowed_cpus=sorted(os.sched_getaffinity(0)),purpose='One actual initial compiled-stage gravity input, reused without another initial density/velocity/gravity fit; no full face-node dump',array_lifetime='one necessary shared mmap input; no per-cell point/design arrays retained in parent; per-worker LRU capacity1')
assert expected_bytes+preflight['other_output_estimate_bytes']<preflight['free_filesystem_bytes'] and expected_bytes+64*1024**3<memory['MemAvailable']/2 and len(preflight['allowed_cpus'])>=4
(args.output/'preflight.json').write_text(json.dumps(preflight,indent=2));print(json.dumps({k:v for k,v in preflight.items() if k!='allowed_cpus'}),flush=True)
node_path=args.output/'gravity_nodes.npy';parent_nodes=np.lib.format.open_memmap(node_path,mode='w+',dtype=float,shape=(ng,7));stage_key=hashlib.sha256(initial_path.read_bytes()).hexdigest()+':'+fitmanifest['source_sha256'][str(native_path)]

def init():
 global nodes,cache
 resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));nodes=np.load(node_path,mmap_mode='r+');cache=VolumeCache(1);cache.begin_stage(stage_key)

def worker(i):
 try:
  row=fit_rows[i];r=raw[i];cell=DiskCell.from_native(cells[i],R);d=CompactPotentialDensity(cell,r,row['prior'],row['parameter'],potential,cache,stage_key,int(ids[i]),24)
  begin=time.perf_counter();s=dict(states[i]);s['density']=d;metadata,diag=prepare_cell(i,s,d,mass[i],factors[i],centers,means,second,tree);s['velocity_quadratic']=metadata;prepare_time=time.perf_counter()-begin
  p=d.points;dm=d.weighted_mass;assert len(p)==offsets[i+1]-offsets[i]
  v=s['v']+(p-s['c'])@s['B'].T+quad.correction(s,p);a=-gravity_coefficient(np.linalg.norm(p,axis=1))[:,None]*p
  begin=time.perf_counter();gravity=source(p,dm,v,a);samples,_=quad.velocity_samples(s,factors[i],mass[i],p,ALPHA)
  pp=dm@v;dd=p.T@(dm[:,None]*v);pi=np.einsum('n,nki,nkj->ij',dm/4,samples,samples);observed=np.r_[pp,dd.ravel(),pi.ravel()];target=np.r_[r[7:9],r[9:13],r[13:17]];defect=observed-target
  if np.any(abs(defect)>1e-14+1e-10*abs(target)):raise RuntimeError(dict(status='prepared_velocity_moment_mismatch',maximum=float(abs(defect).max())))
  gn=np.column_stack([p,dm,v,a]);nodes[offsets[i]:offsets[i+1]]=gn
  gen=model.velocity(s,positions[i]);gravity_time=time.perf_counter()-begin
  # Component absolute sums permit per-cell compiled gravity rounding checks.
  aa=abs(a);xx=abs(p);vv=abs(v);ga=np.zeros(18);ga[7:9]=dm@aa;ga[9:13]=np.einsum('n,ni,nj->ij',dm,xx,aa).ravel();av=np.einsum('n,ni,nj->ij',dm,vv,aa);ga[13:17]=(av+av.T).ravel();ga[17]=np.sum(dm*np.sum(vv*aa,axis=1))
  return dict(id=int(ids[i]),status='passed',metadata={k:(float(x) if np.ndim(x)==0 else np.asarray(x).tolist()) for k,x in metadata.items()},diagnostics=diag,generator=gen.tolist(),gravity=gravity.tolist(),gravity_absolute_sum=ga.tolist(),moment_defect=defect.tolist(),points=len(p),prepare_seconds=prepare_time,gravity_and_moment_seconds=gravity_time,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 except Exception as exc:return dict(id=int(ids[i]),status='failed',error=repr(exc),traceback=traceback.format_exc())
rows=[]
try:
 with mp.get_context('fork').Pool(4,initializer=init) as pool:
  for row in pool.imap(worker,range(N),chunksize=4):
   rows.append(row)
   with (args.output/'cells.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
   if row['status']!='passed':raise RuntimeError(row)
   if len(rows)%500==0:print(json.dumps(dict(cells=len(rows),seconds=time.perf_counter()-started)),flush=True)
 parent_nodes.flush()
 metadata={key:np.array([row['metadata'][key] for row in rows]) for key in rows[0]['metadata']}
 np.savez(args.output/'prepared.npz',ids=ids,state=initial['state'],factors=factors,mass=mass,Psi=initial['Psi'],generator=np.array([r['generator'] for r in rows]),gravity=np.array([r['gravity'] for r in rows]),gravity_absolute_sum=np.array([r['gravity_absolute_sum'] for r in rows]),moment_defect=np.array([r['moment_defect'] for r in rows]),gravity_offsets=offsets,**metadata)
 h=hashlib.sha256()
 with node_path.open('rb') as f:
  while chunk:=f.read(64*1024**2):h.update(chunk)
 report=dict(status='passed_initial_velocity_gravity_preparation',scope=__doc__,cells=N,gravity_nodes=ng,gravity_input_bytes=node_path.stat().st_size,gravity_input_sha256=h.hexdigest(),provider_calls=0,evolution_steps=0,maximum_velocity_moment_defect=max(float(abs(np.array(r['moment_defect'])).max()) for r in rows),maximum_covariance_split_error=max(r['diagnostics']['covariance_split_error'] for r in rows),maximum_orthogonality_error=max(r['diagnostics']['orthogonality'] for r in rows),minimum_velocity_neighbors=min(r['diagnostics']['neighbors'] for r in rows),worker_prepare_wall_sum=sum(r['prepare_seconds'] for r in rows),worker_gravity_and_moment_wall_sum=sum(r['gravity_and_moment_seconds'] for r in rows),peak_worker_RSS_KiB=max(r['peak_RSS_KiB'] for r in rows),seconds=time.perf_counter()-started,initial_sha256=hashlib.sha256(initial_path.read_bytes()).hexdigest(),source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),modules/'compact_potential_density.py',modules/'quadratic_cell_prepare.py',fitroot/'manifest.json',fitroot/'report.json']},limits='Actual initial preparation only; no full closed-root/RHS, dynamic wall flux, timestep evolution or long-time stability claim. Same-state metadata and necessary gravity input are retained for first stage0 reuse.')
except Exception as exc:
 report=dict(status='failed_initial_velocity_gravity_preparation',error=repr(exc),last_row=rows[-1] if rows else None,seconds=time.perf_counter()-started,provider_calls=0,evolution_steps=0)
 (args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(status=report['status'],error=report['error'])),flush=True);raise
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

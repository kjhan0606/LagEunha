"""Bounded CPU preparation with anonymous shared gravity nodes consumed directly.

The per-cell numerical worker is copied from the verified saved-stage worker;
no gravity-node file is written, hashed, reloaded, or retained across stages.
"""
import argparse,hashlib,importlib.util,json,multiprocessing as mp,os,resource,sys,time,traceback
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.spatial import cKDTree
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_full_disk_arc_roots_modules_v7_20261019'
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

import mmap
_context=None

def load_potential():
 spec=importlib.util.spec_from_file_location('in_memory_physical_potential',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 return module.potential
potential=load_potential()

def init():
 global cache
 resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
 cache=VolumeCache(1);cache.begin_stage(_context.stage_key)

def worker(i):
 fit_rows,raw,cells,R,mass,factors,centers,means,second,tree,offsets,ids,positions,states,stage_key=(_context.fit_rows,_context.raw,_context.cells,_context.R,_context.mass,_context.factors,_context.centers,_context.means,_context.second,_context.tree,_context.offsets,_context.ids,_context.positions,_context.states,_context.stage_key)
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

def prepare(initial,native,mesh,fit_rows,stage_key,workers=8):
 global _context,nodes
 started=time.perf_counter()
 assert workers in (4,8) and isinstance(stage_key,str) and stage_key
 assert os.environ['OPENBLAS_NUM_THREADS']==os.environ['OMP_NUM_THREADS']=='1' and os.environ['OMP_PROC_BIND']=='false'
 mass=initial['mass'];ids=initial['ids'];factors=initial['factors']
 raw=np.column_stack([mass,initial['state'][:-4].reshape(-1,20)[:,3:]])
 positions=initial['state'][:-4].reshape(-1,20)[:,:2];N=len(ids)
 assert [r['id'] for r in fit_rows]==ids.tolist() and all(r['status']=='passed' for r in fit_rows)
 assert np.array_equal(mesh['ids'],ids) and np.array_equal(mesh['sites'],positions)
 assert np.array_equal(mesh['weights'],initial['state'][:-4].reshape(-1,20)[:,2])
 model.GAMMA=float(np.float32(1.666667))
 states=[model.decode(None,r,density=SimpleNamespace(mass=float(r[0]),area=float(row['parameter']['area']),origin=np.array(row['parameter']['origin']))) for r,row in zip(raw,fit_rows)]
 for state,L in zip(states,factors):state['R']=L@L.T
 centers=np.array([state['c'] for state in states]);means=np.array([state['v'] for state in states])
 second=raw[:,3:7].reshape(N,2,2)/mass[:,None,None];tree=cKDTree(centers)
 offsets=np.r_[0,np.cumsum([r['nodes24'] for r in fit_rows])].astype(np.int64);ng=int(offsets[-1])
 expected_bytes=ng*7*8
 available=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')))*1024
 assert expected_bytes+workers*16*1024**3<available/2 and len(os.sched_getaffinity(0))>=workers
 # Anonymous MAP_SHARED pages are inherited by workers; each writes one
 # disjoint cell interval. Buffer stays alive through the actual native stage.
 owner=mmap.mmap(-1,expected_bytes,flags=mmap.MAP_SHARED,prot=mmap.PROT_READ|mmap.PROT_WRITE)
 nodes=np.ndarray((ng,7),dtype=float,buffer=owner)
 _context=SimpleNamespace(fit_rows=fit_rows,raw=raw,cells=native['cells'],R=2.5,mass=mass,factors=factors,
  centers=centers,means=means,second=second,tree=tree,offsets=offsets,ids=ids,positions=positions,states=states,stage_key=stage_key)
 rows=[]
 try:
  with mp.get_context('fork').Pool(workers,initializer=init) as pool:
   for row in pool.imap(worker,range(N),chunksize=4):
    if row['status']!='passed':raise RuntimeError(row)
    rows.append(row)
  metadata={key:np.array([row['metadata'][key] for row in rows]) for key in rows[0]['metadata']}
  prepared=dict(ids=ids,state=initial['state'],factors=factors,mass=mass,Psi=initial['Psi'],
   generator=np.array([r['generator'] for r in rows]),gravity=np.array([r['gravity'] for r in rows]),
   gravity_absolute_sum=np.array([r['gravity_absolute_sum'] for r in rows]),
   moment_defect=np.array([r['moment_defect'] for r in rows]),gravity_offsets=offsets,**metadata)
  report=dict(status='passed_in_memory_velocity_gravity_preparation',cells=N,workers=workers,
   gravity_nodes=ng,gravity_input_bytes=expected_bytes,gravity_file_bytes_written=0,
   maximum_velocity_moment_defect=max(float(abs(np.array(r['moment_defect'])).max()) for r in rows),
   maximum_covariance_split_error=max(r['diagnostics']['covariance_split_error'] for r in rows),
   maximum_orthogonality_error=max(r['diagnostics']['orthogonality'] for r in rows),
   seconds=time.perf_counter()-started,peak_worker_RSS_KiB=max(r['peak_RSS_KiB'] for r in rows),
   available_memory_bytes=available,stage_key=stage_key,provider_calls=0,evolution_steps=0)
  return SimpleNamespace(prepared=prepared,gravity_nodes=nodes,owner=owner,report=report)
 finally:
  _context=None
  nodes=None

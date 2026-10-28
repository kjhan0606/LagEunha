"""Same constrained density fit in bounded workers, q48 at checkpoints only.

No fit/root tolerance or q24 rule is changed. Returns compact per-cell fits
in memory and never supplies an analytic profile to an evolved state.
"""
import argparse,hashlib,importlib.util,json,multiprocessing as mp,os,resource,sys,time,traceback
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_full_disk_arc_roots_modules_v7_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_moment_flux_batch_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from curved_disk_cell import DiskCell
sys.path.insert(0,str(base/'Codex_full_disk_fast_modules_v1_20261019'))
from curved_contact_priors_fast import geometric_density_bounds,prior_for_cell
from curved_contact_density import CurvedContactPriorDensity
from contact_potential_density import ContactPotentialDensity,ReconstructionFailure,features
from saved_potential_taylor import taylor_coefficients
from audit_canonical_pair_history import gravity_coefficient

_context=None
ld=np.longdouble
spec=importlib.util.spec_from_file_location('memory_density_potential',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
_potential_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(_potential_module)
potential=_potential_module.potential

def init():resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))

def bounds_worker(i):
 mesh,cells,R,areas,mass,ids=(_context.mesh,_context.cells,_context.R,_context.areas,_context.mass,_context.ids)
 try:
  if mesh['arc_counts'][i]:cell=DiskCell.from_native(cells[i],R)
  else:cell=SimpleNamespace(arcs=[],vertices=[tuple(Fraction(v) for v in edge['A']) for edge in cells[i]['edges']])
  return dict(status='passed',bounds=geometric_density_bounds(cell,float(areas[i]),float(mass[i])))
 except Exception as exc:return dict(status='failed',id=int(ids[i]),error=repr(exc),traceback=traceback.format_exc())

def worker(i):
 cells,R,raw,areas,centers,bounds,neighbors,means,Psi,ids=(_context.cells,_context.R,_context.raw,_context.areas,_context.centers,_context.bounds,_context.neighbors,_context.means,_context.Psi,_context.ids)
 meta=None;parameter=None;timers={}
 try:
  begin=time.perf_counter();cell=DiskCell.from_native(cells[i],R);meta=prior_for_cell(i,cell,raw,areas,centers,bounds,neighbors,means);timers['geometry_and_prior']=time.perf_counter()-begin
  r=raw[i];begin=time.perf_counter();old=CurvedContactPriorDensity(cell,r[0],r[1:3],r[3:7].reshape(2,2),meta,12,(.25,.5));timers['baseline_fit']=time.perf_counter()-begin
  t=taylor_coefficients(old.origin[None,:],potential,gravity_coefficient);s=ld(old.scale);H=t['H'][0].astype(ld)
  reference=np.r_[ld(t['Phi'][0]),s*t['gradient'][0].astype(ld),.5*s*s*H[0,0],s*s*H[0,1],.5*s*s*H[1,1]]
  begin=time.perf_counter();fit=ContactPotentialDensity(cell,r[0],r[1:3],r[3:7].reshape(2,2),meta,old,potential,Psi[i],reference,24,(.25,.5));timers['potential_fit']=time.perf_counter()-begin
  parameter=dict(id=int(ids[i]),theta=fit.theta.tolist(),logZ=fit.logZ,origin=fit.origin.tolist(),scale=fit.scale,area=fit.area,reference=[str(v) for v in fit.reference],projection_remainder=[str(v) for v in fit.remainder_projection],sigma=str(fit.sigma),rank_error=str(fit.rank_error),target=str(fit.Psi_target),residual=fit.residual,iterations=fit.iterations,condition=fit.condition,entropy=fit.entropy)
  if _context.independent_q48:
   begin=time.perf_counter();p,w,prior=fit.quadrature(48);dm=(w*fit._value_with_prior(p,prior)).astype(ld);A=features(p,fit.origin,fit.scale);observed=dm@A;target=ld(r[0])*fit.target_A
   observed_Psi=dm@potential(p.astype(ld));chi=potential(p.astype(ld))-A@fit.reference-A@fit.remainder_projection
   scaled_remainder=(dm@chi-(fit.Psi_target-ld(r[0])*(fit.target_A@fit.projection)))/(ld(r[0])*fit.sigma)
   moment_error=observed-target;Psi_error=observed_Psi-Psi[i]
   gate=bool(np.all(abs(moment_error)<=1e-14+1e-10*abs(target))) and bool(abs(Psi_error)<=1e-14+1e-10*abs(Psi[i]))
   timers['independent_q48']=time.perf_counter()-begin
   return dict(id=int(ids[i]),status='passed' if gate else 'independent_q48_gate_failed',prior=meta,parameter=parameter,normalized_mqQ_error=float(abs(moment_error/ld(r[0])).max()),Psi_error=float(Psi_error),scaled_remainder_error=float(scaled_remainder),nodes24=len(fit.points),nodes48=len(p),minimum_sampled_density=float(dm.min()),timers=timers,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
  timers['independent_q48']=0.
  return dict(id=int(ids[i]),status='passed',prior=meta,parameter=parameter,
   normalized_mqQ_error=None,Psi_error=None,scaled_remainder_error=None,nodes24=len(fit.points),nodes48=0,
   minimum_sampled_density=None,timers=timers,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 except Exception as exc:
  row=dict(id=int(ids[i]),status='failed',error=repr(exc),traceback=traceback.format_exc(),prior=meta,parameter=parameter,timers=timers)
  if isinstance(exc,ReconstructionFailure):row.update(reconstruction_status=exc.status,details=exc.details)
  return row


def prepare(initial,native,mesh,stage_key,workers=8,independent_q48=False):
 global _context
 started=time.perf_counter();assert workers in (4,8) and stage_key
 assert os.environ['OPENBLAS_NUM_THREADS']==os.environ['OMP_NUM_THREADS']=='1' and os.environ['OMP_PROC_BIND']=='false'
 mass,ids,Psi=initial['mass'],initial['ids'],initial['Psi'];N=len(ids)
 assert not initial['frozen'].any() and initial['psi_constrained_mask'].all()
 assert np.array_equal(mesh['ids'],ids)
 aa=initial['state'][:-4].reshape(N,20)
 assert np.array_equal(mesh['sites'],aa[:,:2]) and np.array_equal(mesh['weights'],aa[:,2])
 raw=np.column_stack([mass,aa[:,3:]])
 areas=mesh['moments'][:,0];centers=mesh['moments'][:,1:3]/areas[:,None];means=mass/areas
 lookup={int(g):i for i,g in enumerate(ids)};neighbors=[set() for _ in ids]
 for a,b in mesh['face_keys']:
  i,j=lookup[int(a)],lookup[int(b)];neighbors[i].add(j);neighbors[j].add(i)
 available=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')))*1024
 assert len(os.sched_getaffinity(0))>=workers and workers*16*1024**3<available/2
 _context=SimpleNamespace(mesh=mesh,cells=native['cells'],R=float(mesh['radius']),areas=areas,mass=mass,ids=ids,
  raw=raw,centers=centers,means=means,bounds=None,neighbors=neighbors,Psi=Psi,independent_q48=independent_q48)
 try:
  bounds=[]
  with mp.get_context('fork').Pool(workers,initializer=init) as pool:
   for row in pool.imap(bounds_worker,range(N),chunksize=16):
    if row['status']!='passed':raise RuntimeError(row)
    bounds.append(row['bounds'])
  _context.bounds=bounds;rows=[]
  with mp.get_context('fork').Pool(workers,initializer=init) as pool:
   for row in pool.imap(worker,range(N),chunksize=4):
    if row['status']!='passed':raise RuntimeError(row)
    rows.append(row)
  report=dict(status='passed_in_memory_Psi_density_preparation',cells=N,workers=workers,
   independent_q48=independent_q48,nonflat=sum(not row['prior']['flat'] for row in rows),
   maximum_dual_residual=max(row['parameter']['residual'] for row in rows),
   maximum_condition=max(row['parameter']['condition'] for row in rows),
   minimum_sigma=min(float(row['parameter']['sigma']) for row in rows),
   nodes24=sum(row['nodes24'] for row in rows),nodes48=sum(row['nodes48'] for row in rows),
   worker_wall_timing_sums={key:sum(row['timers'][key] for row in rows) for key in rows[0]['timers']},
   seconds=time.perf_counter()-started,stage_key=stage_key,provider_calls=0,evolution_steps=0)
  return SimpleNamespace(rows=rows,report=report)
 finally:_context=None

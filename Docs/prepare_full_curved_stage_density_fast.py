"""Prepare Psi density on one actual moving full-domain stage.

Same immutable IC/geometry/targets will be consumed by the first actual
provider. No RHS, massroot, evolution, analytical profile injection or fallback.
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
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--state',type=Path,required=True);parser.add_argument('--geometry',type=Path,required=True);parser.add_argument('--workers',type=int,choices=(4,8),default=4);args=parser.parse_args();args.output.mkdir(exist_ok=False);started=time.perf_counter()
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
potential_path=base/'Codex_weight_remap_moments_centered_20261017/audit.py';potential=load('physical_potential',potential_path).potential
initial_path=args.state;initial=np.load(initial_path);mass=initial['mass'];ids=initial['ids'];Psi=initial['Psi'];raw=np.column_stack([mass,initial['state'][:-4].reshape(-1,20)[:,3:]]);N=len(ids)
mesh_path=args.geometry/'mesh.npz';mesh=np.load(mesh_path);assert np.array_equal(mesh['ids'],ids)
geometry_path=args.geometry/'full_exact_supports.json';native=json.loads(geometry_path.read_text());cells=native['cells'];R=float(mesh['radius'])
assert not initial['frozen'].any() and initial['psi_constrained_mask'].all()
areas=mesh['moments'][:,0];centers=mesh['moments'][:,1:3]/areas[:,None];means=mass/areas
lookup={int(g):i for i,g in enumerate(ids)};neighbors=[set() for _ in ids]
for a,b in mesh['face_keys']:
 i=lookup[int(a)];j=lookup[int(b)];neighbors[i].add(j);neighbors[j].add(i)
assert os.environ['OPENBLAS_NUM_THREADS']==os.environ['OMP_NUM_THREADS']=='1' and os.environ['OMP_PROC_BIND']=='false'
ld=np.longdouble
available=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')))*1024
assert len(os.sched_getaffinity(0))>=args.workers and args.workers*16*1024**3<available/2
manifest=dict(scope=__doc__,cells=N,workers=args.workers,available_memory_bytes=available,worker_virtual_memory_limit_bytes=16*1024**3,provider_calls=0,evolution_steps=0,fit_order=24,baseline_preconditioner_order=12,independent_order=48,Psi_mode='transported stage sidecar; all cells constrained',source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),initial_path,mesh_path,geometry_path,potential_path,modules/'contact_potential_density.py',modules/'curved_contact_density.py',modules/'curved_contact_priors.py',modules/'curved_disk_cell.py',base/'Codex_full_disk_fast_modules_v1_20261019/curved_contact_priors_fast.py',base/'Codex_full_disk_fast_modules_v1_20261019/native_exact_prior_bisection.py',base/'Codex_full_disk_fast_modules_v1_20261019/native_exact_prior_bisection.cpp',base/'Codex_full_disk_fast_modules_v1_20261019/libexact_prior.so']},memory_policy='explicitly bounded independent cell workers; one cell and fit/independent quadrature at a time; compact fit metadata only')
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2));np.savez(args.output/'targets.npz',ids=ids,raw=raw,Psi=Psi)

def init():resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
def bounds_worker(i):
 try:
  if mesh['arc_counts'][i]:cell=DiskCell.from_native(cells[i],R)
  else:cell=SimpleNamespace(arcs=[],vertices=[tuple(Fraction(v) for v in edge['A']) for edge in cells[i]['edges']])
  return dict(status='passed',bounds=geometric_density_bounds(cell,float(areas[i]),float(mass[i])))
 except Exception as exc:return dict(status='failed',id=int(ids[i]),error=repr(exc),traceback=traceback.format_exc())

bounds=[]
with mp.get_context('fork').Pool(args.workers,initializer=init) as pool:
 for i,row in enumerate(pool.imap(bounds_worker,range(N),chunksize=16)):
  if row['status']!='passed':
   (args.output/'failure_bounds.json').write_text(json.dumps(row,indent=2));raise RuntimeError(row)
  bounds.append(row['bounds'])
(args.output/'density_bounds.json').write_text(json.dumps(bounds))
print(json.dumps(dict(status='bounds_prepared',cells=N,seconds=time.perf_counter()-started)),flush=True)

def worker(i):
 meta=None;parameter=None;timers={}
 try:
  begin=time.perf_counter();cell=DiskCell.from_native(cells[i],R);meta=prior_for_cell(i,cell,raw,areas,centers,bounds,neighbors,means);timers['geometry_and_prior']=time.perf_counter()-begin
  r=raw[i];begin=time.perf_counter();old=CurvedContactPriorDensity(cell,r[0],r[1:3],r[3:7].reshape(2,2),meta,12,(.25,.5));timers['baseline_fit']=time.perf_counter()-begin
  t=taylor_coefficients(old.origin[None,:],potential,gravity_coefficient);s=ld(old.scale);H=t['H'][0].astype(ld)
  reference=np.r_[ld(t['Phi'][0]),s*t['gradient'][0].astype(ld),.5*s*s*H[0,0],s*s*H[0,1],.5*s*s*H[1,1]]
  begin=time.perf_counter();fit=ContactPotentialDensity(cell,r[0],r[1:3],r[3:7].reshape(2,2),meta,old,potential,Psi[i],reference,24,(.25,.5));timers['potential_fit']=time.perf_counter()-begin
  parameter=dict(id=int(ids[i]),theta=fit.theta.tolist(),logZ=fit.logZ,origin=fit.origin.tolist(),scale=fit.scale,area=fit.area,reference=[str(v) for v in fit.reference],projection_remainder=[str(v) for v in fit.remainder_projection],sigma=str(fit.sigma),rank_error=str(fit.rank_error),target=str(fit.Psi_target),residual=fit.residual,iterations=fit.iterations,condition=fit.condition,entropy=fit.entropy)
  begin=time.perf_counter();p,w,prior=fit.quadrature(48);dm=(w*fit._value_with_prior(p,prior)).astype(ld);A=features(p,fit.origin,fit.scale);observed=dm@A;target=ld(r[0])*fit.target_A
  observed_Psi=dm@potential(p.astype(ld));chi=potential(p.astype(ld))-A@fit.reference-A@fit.remainder_projection
  scaled_remainder=(dm@chi-(fit.Psi_target-ld(r[0])*(fit.target_A@fit.projection)))/(ld(r[0])*fit.sigma)
  moment_error=observed-target;Psi_error=observed_Psi-Psi[i]
  gate=bool(np.all(abs(moment_error)<=1e-14+1e-10*abs(target))) and bool(abs(Psi_error)<=1e-14+1e-10*abs(Psi[i]))
  timers['independent_q48']=time.perf_counter()-begin
  return dict(id=int(ids[i]),status='passed' if gate else 'independent_q48_gate_failed',prior=meta,parameter=parameter,normalized_mqQ_error=float(abs(moment_error/ld(r[0])).max()),Psi_error=float(Psi_error),scaled_remainder_error=float(scaled_remainder),nodes24=len(fit.points),nodes48=len(p),minimum_sampled_density=float(dm.min()),timers=timers,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 except Exception as exc:
  row=dict(id=int(ids[i]),status='failed',error=repr(exc),traceback=traceback.format_exc(),prior=meta,parameter=parameter,timers=timers)
  if isinstance(exc,ReconstructionFailure):row.update(reconstruction_status=exc.status,details=exc.details)
  return row

rows=[];counts={}
try:
 with mp.get_context('fork').Pool(args.workers,initializer=init) as pool:
  for row in pool.imap(worker,range(N),chunksize=4):
   rows.append(row);counts[row['status']]=counts.get(row['status'],0)+1
   with (args.output/'cells.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
   if row['status']!='passed':raise RuntimeError(dict(status='initial_density_failed',id=row['id'],reason=row.get('reconstruction_status',row['error'] if 'error' in row else row['status'])))
   if len(rows)%500==0:print(json.dumps(dict(cells=len(rows),counts=counts,seconds=time.perf_counter()-started)),flush=True)
 report=dict(status='passed_full_stage_Psi_density_preparation',cells=N,provider_calls=0,evolution_steps=0,counts=counts,nonflat=sum(not r['prior']['flat'] for r in rows),maximum_dual_residual=max(r['parameter']['residual'] for r in rows),maximum_condition=max(r['parameter']['condition'] for r in rows),minimum_sigma=min(float(r['parameter']['sigma']) for r in rows),maximum_rank_error=max(float(r['parameter']['rank_error']) for r in rows),maximum_q48_normalized_mqQ_error=max(r['normalized_mqQ_error'] for r in rows),maximum_q48_Psi_error=max(abs(r['Psi_error']) for r in rows),signed_q48_Psi_error=sum(r['Psi_error'] for r in rows),absolute_q48_Psi_error=sum(abs(r['Psi_error']) for r in rows),maximum_q48_scaled_remainder_error=max(abs(r['scaled_remainder_error']) for r in rows),nodes24=sum(r['nodes24'] for r in rows),nodes48=sum(r['nodes48'] for r in rows),worker_wall_timing_sums={name:sum(r['timers'][name] for r in rows) for name in ('geometry_and_prior','baseline_fit','potential_fit','independent_q48')},peak_worker_RSS_KiB=max(r['peak_RSS_KiB'] for r in rows),seconds=time.perf_counter()-started,state_sha256=manifest['source_sha256'][str(initial_path)],time=float(initial['time']),limits='Actual moving-stage fits and independent reconstruction only. No full stage0/root, dynamic boundary ledger, evolution, timestep comparison or long-time success claim. Compact fit parameters retained for reuse, not reconstructed or refitted automatically.')
except Exception as exc:
 report=dict(status='failed_stage_density_preparation',cells_completed=len(rows),counts=counts,error=repr(exc),last_row=rows[-1] if rows else None,seconds=time.perf_counter()-started,provider_calls=0,evolution_steps=0)
 (args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(status=report['status'],error=report['error'])),flush=True);raise
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

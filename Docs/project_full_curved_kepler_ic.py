"""Conservative raw18/Psi/positive residual factors on actual disk cells.

New physical boundary benchmark. Analytical field is used for IC/reference
only. No density fit, fluid provider, massroot, reset, or evolution occurs.
"""
import argparse,hashlib,importlib.util,json,multiprocessing as mp,os,resource,sys,time,traceback
from pathlib import Path
import numpy as np
from scipy.integrate import quad
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_circular_curved_modules_v3_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_steady_ale_reference_modules_20261018'),str(base)]
from curved_disk_cell import DiskCell
from steady_kepler_ale_reference import field_moments
from audit_canonical_pair_history import gravity_coefficient
from run_canonical_av_branches_20261018 import read
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False);started=time.perf_counter()
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
potential_path=base/'Codex_weight_remap_moments_centered_20261017/audit.py';potential=load('actual_softened_potential',potential_path).potential
native_path=base/'Codex_circular_sparse_support_full_20261019/full_exact_supports.json';native=json.loads(native_path.read_text());cells=native['cells'];N=len(cells)
geometry_root=base/'Codex_full_curved_disk_mesh_20261019';geometry_report=json.loads((geometry_root/'report.json').read_text());assert geometry_report['status']=='passed_full_curved_geometry_not_fluid_IC'
mesh=np.load(geometry_root/'mesh.npz');ids=mesh['ids'];sites=mesh['sites'];R=float(mesh['radius']);assert N==len(ids)==19945
assert os.environ['OPENBLAS_NUM_THREADS']==os.environ['OMP_NUM_THREADS']=='1' and os.environ['OMP_PROC_BIND']=='false'
ic,t=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t==0.;P0=float(ic['pressure'][0]);gamma=float(np.float32(1.666667))
env_path=base/'Codex_canonical_initial_20261017/run/environment.json';env=json.loads(env_path.read_text());rin,rout,floor=[float(env['EUNHA_KEPLER_'+k]) for k in ('RIN','ROUT','FLOOR')]
assert rin==.5 and rout==2. and floor==.01 and R==2.5
ld=np.longdouble;eps=np.finfo(float).eps
memory={}
for line in Path('/proc/meminfo').read_text().splitlines():
 k,*v=line.split()
 if k in ('MemTotal:','MemAvailable:'):memory[k[:-1]]=int(v[0])*1024
preflight=dict(cells=N,quad24_nodes=int(mesh['quad24_counts'].sum()),quad48_node_estimate=4*int(mesh['quad24_counts'].sum()),maximum_cell_quad24_nodes=int(mesh['quad24_counts'].max()),workers=4,worker_threads=1,allowed_cpus=sorted(os.sched_getaffinity(0)),worker_virtual_memory_limit_bytes=16*1024**3,memory=memory,estimated_additional_output_bytes=N*3000+10*1024**2,free_output_filesystem_bytes=os.statvfs(base).f_bavail*os.statvfs(base).f_frsize,volume_storage='per worker one cell/order; nodes discarded before next cell')
assert len(preflight['allowed_cpus'])>=4 and 64*1024**3<memory['MemAvailable']/2 and preflight['estimated_additional_output_bytes']<preflight['free_output_filesystem_bytes']
(args.output/'preflight.json').write_text(json.dumps(preflight,indent=2));print(json.dumps({k:v for k,v in preflight.items() if k!='allowed_cpus'}),flush=True)

def init():resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
def projected(cell,order):
 points,w=cell.positive_quadrature(order,(.25,.5,2.));field,velocity=field_moments(points,gravity_coefficient,potential,P0,gamma,rin,rout,floor)
 raw=w.astype(ld)@field[:,:18].astype(ld)
 raw[5]=raw[4];raw[15]=raw[14]
 dm=(w*field[:,0]).astype(ld);Psi=dm@potential(points.astype(ld))
 return raw,Psi,points,dm,velocity,len(w)

def worker(i):
 try:
  c=DiskCell.from_native(cells[i],R);r24,Psi24,x,dm,v,n24=projected(c,24)
  m=dm.sum();center=dm@x.astype(ld)/m;mean=dm@v.astype(ld)/m;dx=x.astype(ld)-center;dv=v.astype(ld)-mean;sqrtm=np.sqrt(dm).astype(float)
  X=sqrtm[:,None]*dx.astype(float);Y=sqrtm[:,None]*dv.astype(float)
  coefficient,residuals,rank,singular=np.linalg.lstsq(X,Y,rcond=None)
  if rank!=2 or singular[-1]<=256*eps*singular[0]:raise RuntimeError('unresolved spatial affine fit rank')
  B=coefficient.T;residual=dv-dx@B.astype(ld).T;weighted=sqrtm[:,None]*residual.astype(float)
  _,upper=np.linalg.qr(weighted,mode='reduced');sg=np.where(np.diag(upper)<0,-1.,1.);L=(sg[:,None]*upper).T
  singularR=np.linalg.svd(weighted,compute_uv=False)
  rank_error=256*eps*float(np.sqrt(m))*(float(np.max(abs(v)))+float(np.linalg.norm(B))*float(np.max(abs(dx))))
  if singularR[-1]<=rank_error or np.any(np.diag(L)<=0):raise RuntimeError(dict(status='unresolved_velocity_residual_rank',singular=singularR.tolist(),error=rank_error))
  directR=np.einsum('n,ni,nj->ij',dm,residual,residual);factorR=L@L.T
  raw=np.asarray(r24,dtype=float);mraw=raw[0];q=raw[1:3];Q=raw[3:7].reshape(2,2);p=raw[7:9];D=raw[9:13].reshape(2,2);Pi=raw[13:17].reshape(2,2)
  S=Q-np.outer(q,q)/mraw;C=D-np.outer(q,p)/mraw;Braw=np.linalg.solve(S,C).T;Rraw=Pi-np.outer(p,p)/mraw-Braw@S@Braw.T
  factor_defect=float(np.max(abs(Rraw-factorR)));factor_bound=512*eps*(np.linalg.norm(Pi,np.inf)+np.linalg.norm(np.outer(p,p)/mraw,np.inf)+np.linalg.norm(Braw@S@Braw.T,np.inf))
  if factor_defect>factor_bound:raise RuntimeError(dict(status='raw_factor_consistency_failed',defect=factor_defect,bound=factor_bound))
  thermal=raw[17]-.5*(raw[13]+raw[16]);assert thermal>0 and np.linalg.eigvalsh(S).min()>0
  # Only compact q24 state survives; q48 is an independent IC check.
  del x,dm,v,dx,dv,X,Y,residual,weighted
  r48,Psi48,x,dm,v,n48=projected(c,48);errors=r24-r48
  gate=bool(np.all(abs(errors)<=1e-14+1e-10*abs(r48))) and bool(abs(Psi24-Psi48)<=1e-14+1e-10*abs(Psi48))
  row=dict(status='passed' if gate else 'independent_q48_gate_failed',id=int(ids[i]),index=i,raw=raw.tolist(),factor=L.ravel().tolist(),Psi=str(Psi24),raw_q24_minus_q48=np.asarray(errors,dtype=float).tolist(),Psi_q24_minus_q48=float(Psi24-Psi48),nodes24=n24,nodes48=n48,minimum_R_eigenvalue=float(np.linalg.eigvalsh(factorR).min()),minimum_S_eigenvalue=float(np.linalg.eigvalsh(S).min()),thermal=float(thermal),raw_factor_defect=factor_defect,raw_factor_roundoff_bound=factor_bound,direct_positive_R_vs_factor=float(np.max(abs(directR-factorR))),rank_singular_min=float(singularR[-1]),rank_error=rank_error,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
  return row
 except Exception as exc:return dict(status='failed',id=int(ids[i]),index=i,error=repr(exc),traceback=traceback.format_exc())

rows=[];counts={}
with mp.get_context('fork').Pool(4,initializer=init) as pool:
 for row in pool.imap(worker,range(N),chunksize=8):
  counts[row['status']]=counts.get(row['status'],0)+1;rows.append(row)
  with (args.output/'cells.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
  if len(rows)%500==0:print(json.dumps(dict(cells=len(rows),counts=counts,seconds=time.perf_counter()-started)),flush=True)
failures=[row for row in rows if row['status']!='passed']
if failures:
 report=dict(status='failed_IC_projection',counts=counts,failed_ids=[r['id'] for r in failures],first=failures[0],seconds=time.perf_counter()-started)
 (args.output/'report.json').write_text(json.dumps(report,indent=2));raise RuntimeError(report)
raw=np.array([r['raw'] for r in rows]);factors=np.array([r['factor'] for r in rows]).reshape(N,2,2);Psi=np.array([r['Psi'] for r in rows],dtype=ld);errors=np.array([r['raw_q24_minus_q48'] for r in rows]);Psi_error=np.array([r['Psi_q24_minus_q48'] for r in rows])
mass=raw[:,0];state=np.r_[np.column_stack([sites,np.zeros(N),raw[:,1:]]).ravel(),np.zeros(4)]
np.savez(args.output/'initial.npz',state=state,factors=factors,mass=mass,ids=ids,Psi=Psi,frozen=np.zeros(N,bool),psi_constrained_mask=np.ones(N,bool),weight_rate_seed=np.zeros(N),time=0.)
np.savez(args.output/'projection_errors.npz',ids=ids,raw_q24_minus_q48=errors,Psi_q24_minus_q48=Psi_error,nodes24=np.array([r['nodes24'] for r in rows]),nodes48=np.array([r['nodes48'] for r in rows]))
# Independent 1-D radial integration of the same unmodified physical field.
breaks=(0.,.25,.5,2.,R)
def density(r):return 1. if rin<=r<=rout else floor
def radial(fun):
 value=0.;error=0.
 for a,b in zip(breaks[:-1],breaks[1:]):
  v,e=quad(fun,a,b,epsabs=1e-13,epsrel=1e-13,limit=200);value+=v;error+=e
 return value,error
kin,kinerr=radial(lambda r:np.pi*density(r)*float(gravity_coefficient(np.array([r]))[0])*r**3)
ang,angerr=radial(lambda r:2*np.pi*density(r)*np.sqrt(float(gravity_coefficient(np.array([r]))[0]))*r**3)
phi,phierr=radial(lambda r:2*np.pi*density(r)*float(potential(np.array([[r,0.]],dtype=ld))[0])*r)
exactmass=np.pi*(floor*R**2+(1-floor)*(rout**2-rin**2));exactQ=np.pi/4*(floor*R**4+(1-floor)*(rout**4-rin**4))
totals=raw.astype(ld).sum(0);energy=kin+np.pi*R**2*P0/(gamma-1);L=totals[10]-totals[11]
expected=np.r_[exactmass,0.,0.,exactQ,0.,0.,exactQ,0.,0.,0.,ang/2,-ang/2,0.,kin,0.,0.,kin,energy]
projection_uncertainty=abs(errors).sum(0);sum_scale=abs(raw).sum(0);roundoff=(128+4*N)*eps*np.maximum(sum_scale,1.)
radial_uncertainty=np.zeros(18);radial_uncertainty[10:12]=angerr/2;radial_uncertainty[[13,16,17]]=kinerr
assert np.all(abs(totals-expected)<=roundoff+projection_uncertainty+radial_uncertainty)
assert abs(Psi.sum()-phi)<=float((128+4*N)*eps*np.sum(abs(Psi)))+abs(Psi_error).sum()+phierr
checkpoint=args.output/'initial.npz'
report=dict(status='passed_new_conservative_disk_IC',scope=__doc__,cells=N,workers=4,provider_calls=0,evolution_steps=0,radius=R,physics=dict(rin=rin,rout=rout,floor=floor,P0=P0,gamma=gamma,Re='infinity',physical_viscosity=0.),counts=counts,total_mass=float(totals[0]),exact_mass=exactmass,total_origin_L=float(L),radial_origin_L=ang,total_energy=float(totals[17]),radial_energy=energy,total_Psi=float(Psi.sum()),radial_Psi=phi,total_raw_minus_radial=np.asarray(totals-expected,dtype=float).tolist(),raw_total_roundoff_bound=roundoff.tolist(),raw_q24_minus_q48_absolute_sum=projection_uncertainty.tolist(),Psi_q24_minus_q48_absolute_sum=float(abs(Psi_error).sum()),minimum_factor_R=min(r['minimum_R_eigenvalue'] for r in rows),minimum_thermal=min(r['thermal'] for r in rows),maximum_raw_factor_defect=max(r['raw_factor_defect'] for r in rows),maximum_direct_positive_R_vs_factor=max(r['direct_positive_R_vs_factor'] for r in rows),minimum_rank_margin=min(r['rank_singular_min']/r['rank_error'] for r in rows),quad24_nodes=sum(r['nodes24'] for r in rows),quad48_nodes=sum(r['nodes48'] for r in rows),worker_maximum_RSS_KiB=max(r['peak_RSS_KiB'] for r in rows),seconds=time.perf_counter()-started,checkpoint_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),modules/'curved_disk_cell.py',native_path,geometry_root/'report.json',potential_path,env_path,Path(sys.modules[gravity_coefficient.__module__].__file__),Path(sys.modules[field_moments.__module__].__file__)]},limits='New disk IC and independent quadrature/radial checks only; not old IC bitwise reproduction or self-consistent fluid evolution. Direct positive residual covariance is preserved without floor/noise; raw-factor differences reported. Full-domain density reconstruction, closed root and prescribed two dt evolution remain outstanding.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','raw_total_roundoff_bound','raw_q24_minus_q48_absolute_sum')},indent=2))

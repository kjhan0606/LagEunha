"""Compare the actual one-step HG/RK4/HG map to its response generator."""
import argparse, json, os, re, shutil, subprocess, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read, sha
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');src=b/'Codex_coupled_onset_to_t4_20261018'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root',type=Path,default=b/'Codex_coupled_generator_steps_20261018')
parser.add_argument('--threads',type=int,choices=[1,4],default=4)
parser.add_argument('--dt',type=float,nargs='+',default=[1e-4,5e-5])
parser.add_argument('--dump-stages',action='store_true',help='Dump first six force calls of the single split step.')
parser.add_argument('--binary',type=Path,default=b/'Codex_orbit_work_20261018/eunha2.exe')
parser.add_argument('--smooth-qr',action='store_true')
parser.add_argument('--adaptive-qr',action='store_true')
parser.add_argument('--reject-once',action='store_true',help='Reject a trial at twice dt, then accept dt.')
args=parser.parse_args();root=args.root
if args.adaptive_qr:args.smooth_qr=True
assert len(set(args.dt))==len(args.dt) and all(0<dt<=1e-4 for dt in args.dt)
budget=json.loads((b/'Codex_coupled_response_forces_t3p4_20261018/analysis.json').read_text())
initial={};configs={}
for name in ['plus','minus']:
 a,t=read(src/name/'khout.000743.dat');assert abs(t-budget['time'])<1e-12
 initial[name]=a;configs[name]=json.loads((src/name/'environment.json').read_text())
assert shutil.disk_usage(b).free>1024**3
root.mkdir(exist_ok=False);jobs=[];manifest={'source_time':t,'cases':{}}
for dt in args.dt:
 for sign in ['plus','minus']:
  name=f'{sign}_{dt:g}';d=root/name;d.mkdir()
  for f in ['params.dat','khout.000743.dat']:shutil.copyfile(src/sign/f,d/f)
  c=configs[sign].copy();c.update(HYDRO_TSTOP=str(t+dt/2),GFS_MAX_DT=str(dt),EUNHA_DUMP_DT='1000000',EUNHA_DUMP_FINAL='1',OMP_NUM_THREADS=str(args.threads))
  if args.dump_stages:c.update(GFS_TCF_DUMP=str(d/'rhs'),GFS_TCF_DUMP_CALL='1',GFS_TCF_DUMP_LAST='6')
  if args.smooth_qr:c['GFS_QR_SMOOTH']='1'
  if args.adaptive_qr:c['GFS_QR_ADAPTIVE']='1'
  assert not any(int(c.get(k,'0')) for k in ['GFS_MPI_GEOMETRY_ONLY','GFS_TEST_REJECT_FINAL_ONCE'])
  if args.reject_once:c.update(GFS_TEST_REJECT_FINAL_ONCE='1',GFS_MAX_DT=str(2*dt))
  (d/'environment.json').write_text(json.dumps(c,indent=2))
  env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(c)
  f=(d/'log').open('x');p=subprocess.Popen(['timeout','600','mpirun','-np','1',str(args.binary),'params.dat','743'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
  jobs.append((name,d,p,f));manifest['cases'][name]={'dt':dt,'sign':sign,'pid':p.pid,'source_sha256':sha(d/'khout.000743.dat')};print(name,p.pid,flush=True)
(root/'manifest.json').write_text(json.dumps(manifest,indent=2))
codes=[]
for name,d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));codes.append(code);print(name,code,flush=True)
assert all(code==0 for code in codes)
ic,_=read(b/'Codex_canonical_initial_20261017/run/khout.000000.dat');r0=np.hypot(ic['x']-3,ic['y']-3)
masks={'inner':r0<.5,'shear_band':(r0>=.3)&(r0<.5),'central':r0<.2,'edge':(r0>=.5)&(r0<.6),'core':(r0>=.6)&(r0<1.9)}
def energy(states):
 def u(a):
  x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']))
  return v-np.sqrt(gravity_coefficient(np.linalg.norm(x,axis=1)))[:,None]*np.column_stack((-x[:,1],x[:,0]))
 du=(u(states['plus'])-u(states['minus']))/.002
 return {k:float(.5*np.sum(ic['mass'][mask,None]*du[mask]**2)) for k,mask in masks.items()}
e0=energy(initial);result={}
for dt in args.dt:
 end={}
 for sign in ['plus','minus']:
  d=root/f'{sign}_{dt:g}';a,te=read(d/'khout.000744.dat')
  assert abs(te-t-dt)<1e-12
  assert np.array_equal(a['id'],ic['id']) and np.array_equal(a['mass'],ic['mass'])
  rows=[dict(re.findall(r'(\w+)=([^\s]+)',l)) for l in (d/'log').read_text().splitlines() if l.startswith('[ACCEPTED_STATE]')]
  assert len(rows)==1 and rows[0]['finite']=='1' and rows[0]['retries']==str(int(args.reject_once))
  assert abs(float(rows[0]['dt'])-dt)<1e-15
  end[sign]=a
 e1=energy(end);result[str(dt)]={}
 for k in masks:
  derivative=(e1[k]-e0[k])/dt
  result[str(dt)][k]={'measured_forward_derivative':derivative}
  if not args.smooth_qr:
   pred=budget['regions'][k]['total_generator']
   result[str(dt)][k].update(instantaneous_generator=pred,relative_difference=(derivative-pred)/abs(pred))
out={'time':t,'comparisons':result,'scope':'Actual single accepted HG/RK4/HG steps versus zero-step response generator. Forward derivative error is expected O(dt); this is not a formal RK4 convergence test or a long-time accuracy pass.'}
if args.smooth_qr:out['scope']='Smooth QR one-step forward response derivatives only. Legacy generator is intentionally not used for this changed operator; smooth-generator closure and long-time validation remain pending.'
out['binary_sha256']=sha(args.binary);out['smooth_qr']=args.smooth_qr
out['adaptive_qr']=args.adaptive_qr
out['forced_final_rejection']=args.reject_once
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))

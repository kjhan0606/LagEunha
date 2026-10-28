"""Prepare isentropic position/velocity perturbations; no evolution."""
import json,sys,shutil,argparse
from pathlib import Path
import numpy as np
from scipy.spatial import Voronoi
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');src=b/'Codex_onset_velocity_tangent_20261018/reference'
parser=argparse.ArgumentParser()
parser.add_argument('--root',type=Path,default=b/'Codex_coupled_onset_seed_20261018')
parser.add_argument('--epsilon',type=float,default=.001)
args=parser.parse_args();assert 0<args.epsilon<=.001
root=args.root;root.mkdir(exist_ok=False)
layout=json.loads((b/'Codex_pressure_seed_20261018/layout.json').read_text())
fields=['id','x','y','vx','vy','mass','volume','den','pressure','ie','csound','stress.K']
dtype=np.dtype({'names':fields,'formats':['<u8']+['<f8']*(len(fields)-1),'offsets':[layout['u4if']]+[layout[k] for k in fields[1:]],'itemsize':layout['size']})
raw=(src/'khout.000184.dat').read_bytes();a=np.frombuffer(raw,dtype=dtype,count=36864,offset=16)
ids=(a['id']&layout['id_mask']).astype(int);assert np.array_equal(np.sort(ids),np.arange(len(a)))
_,time=read(src/'khout.000184.dat');ic,_=read(b/'Codex_canonical_initial_20261017/run/khout.000000.dat');ic=ic[ids]
x0=np.column_stack((ic['x']-3,ic['y']-3));om=np.sqrt(gravity_coefficient(np.linalg.norm(x0,axis=1)));th=om*time
xc=np.column_stack((x0[:,0]*np.cos(th)-x0[:,1]*np.sin(th),x0[:,0]*np.sin(th)+x0[:,1]*np.cos(th)))
vc=om[:,None]*np.column_stack((-xc[:,1],xc[:,0]));xc+=3
x=np.column_stack((a['x'],a['y']));v=np.column_stack((a['vx'],a['vy']));r=np.linalg.norm(x-3,axis=1)
dx=x-xc;dv=v-vc;mask=r<.6
for seed in [dx,dv]:
 seed[~mask]=0;seed[mask]-=np.average(seed[mask],axis=0,weights=a['mass'][mask])
inside=np.flatnonzero(r<2.4)
def volumes(pos):
 mesh=Voronoi(pos);out=np.empty(len(inside))
 for k,i in enumerate(inside):
  reg=mesh.regions[mesh.point_region[i]];assert reg and min(reg)>=0
  q=mesh.vertices[reg]-pos[i]
  out[k]=abs(np.dot(q[:,0],np.roll(q[:,1],-1))-np.dot(q[:,1],np.roll(q[:,0],-1)))/2
 return out
V0=volumes(x);err=float(np.max(abs(V0/a['volume'][inside]-1)));assert err<1e-9
gamma=float(np.float32(1.666667));cases={}
for name,eps in [('reference',0.),('plus',args.epsilon),('minus',-args.epsilon)]:
 d=root/name;d.mkdir();buf=bytearray(raw);q=np.frombuffer(buf,dtype=dtype,count=len(a),offset=16)
 ratio=np.ones(len(a))
 if eps:
  pos=x+eps*dx;ratio[inside]=volumes(pos)/V0
  q['x']=pos[:,0];q['y']=pos[:,1];q['vx']+=eps*dv[:,0];q['vy']+=eps*dv[:,1]
  q['volume']*=ratio;q['den']/=ratio;q['pressure']*=ratio**(-gamma)
  q['ie']*=ratio**(1-gamma);q['csound']*=ratio**((1-gamma)/2)
 assert np.all(q['volume']>0) and np.all(q['ie']>0)
 assert np.array_equal(q['mass'],a['mass']) and np.array_equal(q['stress.K'],a['stress.K'])
 ent=np.max(abs((q['pressure']/a['pressure'])/(q['den']/a['den'])**gamma-1));assert ent<1e-12
 remote_change=np.max(abs(ratio[r>1]-1));assert remote_change<1e-9
 (d/'khout.000184.dat').write_bytes(buf);shutil.copyfile(src/'params.dat',d/'params.dat')
 c=json.loads((src/'environment.json').read_text());c.update(OMP_NUM_THREADS='4',GFS_ORBIT_WORK_AUDIT='1')
 (d/'environment.json').write_text(json.dumps(c,indent=2))
 cases[name]={'epsilon':eps,'max_relative_volume_change':float(np.max(abs(ratio-1))),'entropy_relative_error':float(ent),'r_gt_1_volume_change':float(remote_change),'snapshot_sha256':sha(d/'khout.000184.dat')}
out={'time':time,'source':str(src/'khout.000184.dat'),'source_sha256':sha(src/'khout.000184.dat'),'baseline_geometry_max_relative_error':err,'cases':cases,'scope':'Initial diagnostic branches only. Position and velocity seeds follow actual departure from exact circular trajectories inside current r<.6, translation removed. Rebuilt volume ratios and fixed cell entropy/mass determine thermodynamics. No entropy perturbation; not a full eigenvector. Solver initialization consistency still requires a check before interpreting evolved response.'}
(root/'manifest.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

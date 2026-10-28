"""Linear pressure-force decomposition at frozen actual checkpoint geometry."""
import sys,json,argparse
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
base=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser()
parser.add_argument('--root',type=Path,default=base/'Codex_qbp_t10_decomposition_20261018')
parser.add_argument('--label',default='t10')
args=parser.parse_args();p=args.root
s=np.loadtxt(p/'raw_faces.dat.rhs.rank0000');s=s[np.argsort(s[:,0])]
f=np.loadtxt(p/'raw_faces.dat.common.rank0000');q=np.loadtxt(p/'raw_faces.dat.rhs.faces.rank0000')
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
assert np.array_equal(s[:,0],np.arange(len(ic))) and np.array_equal(s[:,7],ic['mass'])
g=float(np.float32(1.666667));P=s[:,5];Piso=ic['pressure']*(ic['volume']/s[:,6])**g;Pent=P-Piso
i,j=f[:,:2].astype(int).T;dx=f[:,2:4];A=f[:,4];d=np.linalg.norm(dx,axis=1);e=dx/d[:,None];t=f[:,5:7]-.5*dx
t-=np.sum(t*e,axis=1)[:,None]*e
def force(pr):
 z=-A[:,None]*(.5*(pr[i]+pr[j])[:,None]*e+((pr[j]-pr[i])/d)[:,None]*t)
 return np.column_stack([np.bincount(i,weights=z[:,k],minlength=len(s)) for k in range(2)])
parts={k:force(v) for k,v in [('actual',P),('isentropic',Piso),('entropy_excess',Pent)]}
fq=np.zeros((len(s),2));np.add.at(fq,q[:,0].astype(int),-(q[:,4]*q[:,14])[:,None]*q[:,6:8])
x=s[:,1:3]-3;r=np.linalg.norm(x,axis=1);vc=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack((-x[:,1],x[:,0]));dv=s[:,3:5]-vc
r0=np.hypot(ic['x']-3,ic['y']-3);masks={'inner':(r0<.5)&(ic['den']<.5),'edge':(r0<.6)&(ic['den']>.5),'core':(r0>=.6)&(r0<1.9)&(ic['den']>.5)}
out={'label':args.label,'source_root':str(p),
 'scope':'Frozen evolved state. P_iso=P_initial*(V_initial/V)^gamma at fixed mass; P_entropy=P-P_iso. Linear pressure-force split only, not a dynamically evolved counterfactual or causal fraction.',
 'force_split_max_abs':float(np.max(abs(parts['actual']-parts['isentropic']-parts['entropy_excess']))),'regions':{}}
active=r<1.9
out['pressure_force_reproduction_relative_L2']=float(np.linalg.norm((parts['actual']-(s[:,8:10]-fq))[active])/np.linalg.norm((s[:,8:10]-fq)[active]))
assert out['pressure_force_reproduction_relative_L2']<1e-4
for name,mask in masks.items():
 row={}
 for k,F in parts.items():
  row[k]={'work_on_velocity_departure':float(np.sum(dv[mask]*F[mask])),
    'work_on_circular_velocity':float(np.sum(vc[mask]*F[mask])),
    'acceleration_mass_RMS':float(np.sqrt(np.average(np.sum((F[mask]/s[mask,7,None])**2,axis=1),weights=s[mask,7])))}
 out['regions'][name]=row
with (p/'pressure_entropy_force.json').open('x') as fp:json.dump(out,fp,indent=2)
print(json.dumps(out,indent=2))

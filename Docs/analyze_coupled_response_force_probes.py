"""Instantaneous continuous-generator response budget, including HG modes.

HG modal sums are the zero-step generator, not finite ordered split updates.
"""
import argparse, json, struct, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read, sha
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_coupled_response_forces_t3p4_20261018'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--fields-only',action='store_true')
parser.add_argument('--root',type=Path,default=root)
parser.add_argument('--case-prefix',default='')
parser.add_argument('--output',default='analysis.json')
args=parser.parse_args()
root=args.root
source=b/'Codex_coupled_onset_to_t4_20261018'
kin=json.loads((source/'kinematic_response_000743.json').read_text())
ic,_=read(b/'Codex_canonical_initial_20261017/run/khout.000000.dat')
states=[];hydro=[];hg=[];checks={};eps=kin['epsilon']
for name in ['plus','minus']:
 d=root/(args.case_prefix+name);assert (d/'exit.code').read_text().strip()=='0'
 log=(d/'log').read_text();assert '[MPI_FORCE_ONLY]' in log and 'HG_apply=0' in log
 assert '[ACCEPTED_STATE]' not in log
 if (d/'manifest.json').exists():meta=json.loads((d/'manifest.json').read_text())
 else:
  meta=json.loads((root/'manifest.json').read_text())['cases'][args.case_prefix+name].copy()
  meta['source_copy']='khout.000743.dat'
 p=d/meta['source_copy']
 assert sha(p)==meta['source_sha256']
 assert sha(p)==sha(source/name/'khout.000743.dat'), 'Reused kinematic/pressure budget requires identical states'
 a,t=read(p);assert abs(t-kin['time'])<1e-12
 assert np.array_equal(a['id'],ic['id']) and np.array_equal(a['mass'],ic['mass'])
 states.append(a);n=len(a)
 state=np.loadtxt(d/'raw_faces.dat.hgstate.rank0000');state=state[np.argsort(state[:,0])]
 rhs=np.loadtxt(d/'raw_faces.dat.rhs.rank0000');rhs=rhs[np.argsort(rhs[:,0])]
 assert state.shape==(n,13) and rhs.shape==(n,25)
 assert np.array_equal(state[:,0],np.arange(n)) and np.array_equal(rhs[:,0],np.arange(n))
 assert np.array_equal(a['id'] & 72057594037927935,np.arange(n))
 for field,col in [('x',3),('y',4),('mass',5),('vx',6),('vy',7),('ie',8)]:
  assert np.max(abs(state[:,col]-a[field]))<1e-12
 for field,col in [('x',1),('y',2),('vx',3),('vy',4),('mass',7)]:
  assert np.array_equal(rhs[:,col],a[field])
 interior=np.hypot(a['x']-3,a['y']-3)<1.9
 initerr=max(float(np.max(abs(rhs[interior,col]/a[field][interior]-1))) for field,col in [('volume',6),('pressure',5)])
 assert initerr<1e-9
 hydro.append(rhs[:,8:10])
 path=d/'raw_faces.dat.hgm.rank0000'
 with path.open('rb') as f:magic,nm,ns,ncol=struct.unpack('<4q',f.read(32))
 assert magic==20261017 and ncol==9 and path.stat().st_size==32+nm*72+ns*24
 m=np.memmap(path,mode='r',dtype='<f8',offset=32,shape=(nm,9))
 s=np.memmap(path,mode='r',dtype='<f8',offset=32+nm*72,shape=(ns,3))
 assert np.isfinite(m).all() and np.isfinite(s).all()
 owner=m[:,0].astype(int);start=m[:,2].astype(int);length=m[:,3].astype(int);ids=s[:,0].astype(int)
 assert np.all(m[:,4:6]==0) and np.all(length>0)
 assert np.array_equal(start,np.r_[0,np.cumsum(length[:-1])]) and length.sum()==ns
 assert np.all((ids>=0)&(ids<n)) and np.all((owner>=0)&(owner<n))
 mass=state[:,5];U=state[:,8];v=state[:,6:8]
 frozen=np.linalg.norm(state[:,3:5]-3,axis=1)>=2.5;free=~frozen[ids]
 alpha=np.add.reduceat(np.sum(s[:,1:3]**2,axis=1)/mass[ids]*free,start)
 entry_dot=np.sum(s[:,1:3]*v[ids],axis=1);residual=np.add.reduceat(entry_dot,start)
 C=m[:,6]*(m[:,7]*abs(residual)+np.sqrt((m[:,7]*residual)**2+m[:,8]*U[owner]));C=np.where(alpha>0,C,0.)
 assert np.isfinite(C).all() and np.all(C>=0)
 force=np.zeros((n,2));np.add.at(force,ids,-np.repeat(C*residual,length)[:,None]*s[:,1:3]*free[:,None]);hg.append(force)
 heat=float(np.sum(C[~frozen[owner]]*residual[~frozen[owner]]**2))
 boundary_work=float(np.sum(C*residual*np.add.reduceat(entry_dot*(~free),start)))
 boundary_heat=float(np.sum(C[frozen[owner]]*residual[frozen[owner]]**2))
 kinetic=float(np.sum(v*force));closure=kinetic+heat-boundary_work+boundary_heat
 assert abs(closure)<1e-9*max(heat,abs(kinetic),1e-30)
 checks[name]={'mode_count':int(nm),'HG_energy_closure':closure,'initialization_relative_error':initerr,'snapshot_sha256':sha(p)}
def u(a):
 x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']))
 return v-np.sqrt(gravity_coefficient(np.linalg.norm(x,axis=1)))[:,None]*np.column_stack((-x[:,1],x[:,0]))
du=(u(states[0])-u(states[1]))/(2*eps)
dF=(hydro[0]-hydro[1])/(2*eps);dHG=(hg[0]-hg[1])/(2*eps)
r0=np.hypot(ic['x']-3,ic['y']-3)
masks={'inner':r0<.5,'shear_band':(r0>=.3)&(r0<.5),'central':r0<.2,'edge':(r0>=.5)&(r0<.6),'core':(r0>=.6)&(r0<1.9)}
rows={}
for name,mask in masks.items():
 row=kin['regions'][name].copy();pressure=row['reversible_pressure']
 hydro_power=float(np.sum(du[mask]*dF[mask]));hg_power=float(np.sum(du[mask]*dHG[mask]))
 row.update(hydro_RHS=hydro_power,hydro_minus_reversible_pressure=hydro_power-pressure,HG_generator=hg_power,
            total_generator=hydro_power+hg_power+row['gravity_plus_reference'])
 rows[name]=row
out={'time':kin['time'],'regions':rows,'checks':checks,
 'scope':'Instantaneous finite central-response generator budget using actual production hydro RHS, reconstructed HG modal generator, analytic gravity and circular-reference Jacobian. Hydro-minus-pressure includes TCF and any other hydro RHS terms; no unverified purity claim. HG generator differs from finite split-map derivative. Not an integrated accepted-step budget or convergence/eigenmode proof.'}
if args.fields_only:
 with (root/'force_fields.npz').open('xb') as f:
  np.savez(f,hydro_plus=hydro[0],hydro_minus=hydro[1],hg_plus=hg[0],hg_minus=hg[1])
else:
 with (root/args.output).open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))

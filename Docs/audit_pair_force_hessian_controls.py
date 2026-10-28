"""Force-Jacobian check of thermal Hessian with fixed-entropy controls."""
import sys,json
from pathlib import Path
import numpy as np
from scipy.spatial import Voronoi
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
base=Path('/gpfs/kjhan/LagForce');s=np.loadtxt(base/'Codex_qbp_t10_decomposition_20261018/raw_faces.dat.rhs.rank0000');s=s[np.argsort(s[:,0])]
prior=json.loads((base/'Codex_pair_thermal_curvature_20261018/hessian.json').read_text())
i,j=prior['pair'];x=s[:,1:3];gamma=float(np.float32(1.666667));ids=np.where(np.linalg.norm(x-3,axis=1)<2.4)[0]
normal=(x[j]-x[i])/np.linalg.norm(x[j]-x[i]);basis=np.array([normal,[-normal[1],normal[0]]]);gap=np.linalg.norm(x[j]-x[i])
ci=-s[j,7]/(s[i,7]+s[j,7]);cj=s[i,7]/(s[i,7]+s[j,7])
def geometry(z):
 mesh=Voronoi(z);V=np.full(len(s),np.nan)
 for k in ids:
  reg=mesh.regions[mesh.point_region[k]];assert reg and min(reg)>=0
  q=mesh.vertices[reg]-z[k];V[k]=abs(np.dot(q[:,0],np.roll(q[:,1],-1))-np.dot(q[:,1],np.roll(q[:,0],-1)))/2
 faces=[]
 for pair,vertices in zip(mesh.ridge_points,mesh.ridge_vertices):
  if i not in pair and j not in pair:continue
  assert len(vertices)==2 and min(vertices)>=0
  a,b=pair;ends=mesh.vertices[vertices];A=np.linalg.norm(ends[1]-ends[0]);dx=z[b]-z[a];d=np.linalg.norm(dx);e=dx/d;t=ends.mean(axis=0)-.5*(z[a]+z[b]);t-=np.dot(t,e)*e
  faces.append((a,b,A,d,e,t))
 return V,faces
V0,_=geometry(x);ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
iso=ic['pressure']*(ic['volume']/V0)**gamma
models={'actual':s[:,5],'uniform':np.full(len(s),ic['pressure'][0]),'initial_entropy':iso,'entropy_excess':s[:,5]-iso}
def projected_force(P,faces):
 Fi=np.zeros(2);Fj=np.zeros(2)
 for a,b,A,d,e,t in faces:
  f=-A*(.5*(P[a]+P[b])*e+(P[b]-P[a])*t/d)
  if a==i:Fi+=f
  if b==i:Fi-=f
  if a==j:Fj+=f
  if b==j:Fj-=f
 return basis@(ci*Fi+cj*Fj)
out={'pair':[i,j],'scope':'Restricted pair force derivative, same t10 geometry, entropy fixed separately within each pressure model; no gravity/rotation/dissipation.','samples':[]}
for frac in [1e-3,5e-4]:
 eps=frac*gap;H={k:np.empty((2,2)) for k in models}
 for axis,vec in enumerate(basis):
  motion=np.zeros_like(x);motion[i]=ci*vec;motion[j]=cj*vec
  vp,fp=geometry(x+eps*motion);vm,fm=geometry(x-eps*motion)
  for name,P in models.items():
   qp=projected_force(P*(V0/vp)**gamma,fp);qm=projected_force(P*(V0/vm)**gamma,fm)
   H[name][:,axis]=-(qp-qm)/(2*eps)
 row={'epsilon_over_gap':frac,'models':{}}
 for name,h in H.items():
  row['models'][name]={'matrix':h.tolist(),'antisymmetric_max_abs':float(np.max(abs(h-h.T))), 'symmetric_eigenvalues':np.linalg.eigvalsh(.5*(h+h.T)).tolist()}
 old=np.array(next(v['hessian'] for v in prior['samples'] if v['epsilon_over_gap']==frac))
 row['actual_energy_hessian_max_abs_difference']=float(np.max(abs(H['actual']-old)))
 row['actual_minus_iso_minus_entropy_max_abs']=float(np.max(abs(H['actual']-H['initial_entropy']-H['entropy_excess'])))
 out['samples'].append(row)
root=base/'Codex_pair_force_hessian_controls_20261018';root.mkdir(exist_ok=False)
(root/'result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

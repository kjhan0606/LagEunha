"""Pressure forcing from exact-orbit geometry versus actual early states."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.spatial import Voronoi
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_circular_pressure_forcing_v2_20261018';root.mkdir(exist_ok=False)
ic,_=read(b/'Codex_canonical_initial_20261017/run/khout.000000.dat');n=len(ic);g=float(np.float32(1.666667))
x0=np.column_stack((ic['x']-3,ic['y']-3));r0=np.hypot(ic['x']-3,ic['y']-3);om=np.sqrt(gravity_coefficient(r0));active=r0<2.5
masks={'inner':(r0<.5)&(ic['den']<.5),'edge':(r0<.6)&(ic['den']>.5),'core':(r0>=.6)&(r0<1.9)&(ic['den']>.5)}
def geometry(x):
 mesh=Voronoi(x);V=np.full(n,np.nan)
 for i in range(n):
  reg=mesh.regions[mesh.point_region[i]]
  if not reg or min(reg)<0:continue
  q=mesh.vertices[reg]-x[i];V[i]=abs(np.dot(q[:,0],np.roll(q[:,1],-1))-np.dot(q[:,1],np.roll(q[:,0],-1)))/2
 pairs=[];vertices=[]
 for ij,rr in zip(mesh.ridge_points,mesh.ridge_vertices):
  if len(rr)==2 and min(rr)>=0 and np.isfinite(V[ij]).all():pairs.append(ij);vertices.append(rr)
 i,j=np.asarray(pairs).T;q=mesh.vertices[np.asarray(vertices)];A=np.linalg.norm(q[:,1]-q[:,0],axis=1)
 dx=x[j]-x[i];d=np.linalg.norm(dx,axis=1);e=dx/d[:,None];t=q.mean(axis=1)-.5*(x[i]+x[j]);t-=np.sum(t*e,axis=1)[:,None]*e
 def force(P):
  z=-A[:,None]*(.5*(P[i]+P[j])[:,None]*e+((P[j]-P[i])/d)[:,None]*t)
  F=np.zeros((n,2));np.add.at(F,i,z);np.add.at(F,j,-z);return F
 return V,force
out={'scope':'Independent finite Voronoi faces, compare by persistent ID in lab frame. Reference pressure is isentropic on prescribed exact circular geometry. Frozen-state forcing comparison, not a trajectory counterfactual or additive causal percentage. Only interior material regions are interpreted.','cases':{}}
for label,step in [('t2',184),('t4',368)]:
 a,time=read(b/f'Codex_HG_constrained_20261017/omp8_t20/khout.{step:06d}.dat')
 assert np.array_equal(a['mass'],ic['mass'])
 x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']));r=np.linalg.norm(x,axis=1)
 V,force=geometry(x);sel=r<1.9;verr=float(np.max(abs(V[sel]/a['volume'][sel]-1)));assert verr<1e-9
 Fa=force(a['pressure']);Fi=force(ic['pressure']*(ic['volume']/V)**g)
 theta=om*time;xc=x0.copy();xc[active,0]=(x0[:,0]*np.cos(theta)-x0[:,1]*np.sin(theta))[active];xc[active,1]=(x0[:,0]*np.sin(theta)+x0[:,1]*np.cos(theta))[active]
 Vc,fc=geometry(xc);Pc=ic['pressure']*(ic['volume']/Vc)**g;Fc=fc(Pc)
 du=v-np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack((-x[:,1],x[:,0]))
 row={'time':time,'actual_volume_relative_error':verr,'regions':{}}
 if label=='t2':
  p=b/'Codex_pair_hg_modes_t2_20261018';s=np.loadtxt(p/'raw_faces.dat.rhs.rank0000');s=s[np.argsort(s[:,0])];q=np.loadtxt(p/'raw_faces.dat.rhs.faces.rank0000')
  fq=np.zeros((n,2));np.add.at(fq,q[:,0].astype(int),-(q[:,4]*q[:,14])[:,None]*q[:,6:8]);saved=s[:,8:10]-fq
  err=float(np.linalg.norm((Fa-saved)[sel])/np.linalg.norm(saved[sel]));assert err<1e-4
  row['pressure_force_reproduction_relative_L2']=err
 for name,mask in masks.items():
  vals={}
  for key,F in [('actual',Fa),('actual_isentropic',Fi),('exact_orbit_isentropic',Fc),('actual_minus_reference',Fa-Fc)]:
   vals[key]={'work_on_actual_departure':float(np.sum(du[mask]*F[mask])),
    'acceleration_mass_RMS':float(np.sqrt(np.average(np.sum((F[mask]/a['mass'][mask,None])**2,axis=1),weights=a['mass'][mask])))}
  row['regions'][name]=vals
 out['cases'][label]=row
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))

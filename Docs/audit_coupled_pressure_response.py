"""Pressure contribution to finite-perturbation departure-energy rate."""
import json,sys,argparse
from pathlib import Path
import numpy as np
from scipy.spatial import Voronoi
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
ap=argparse.ArgumentParser();ap.add_argument('--step',type=int,required=True);args=ap.parse_args()
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_coupled_onset_to_t4_20261018';seed=b/'Codex_coupled_onset_seed_20261018'
initial,t0=read(seed/'reference/khout.000184.dat');ic,_=read(b/'Codex_canonical_initial_20261017/run/khout.000000.dat');states=[];eps=.001;gamma=float(np.float32(1.666667));n=len(initial)
for name in ['plus','minus']:
 a,time=read(root/name/f'khout.{args.step:06d}.dat');assert np.array_equal(a['id'],initial['id']);states.append(a)
def geometry(a):
 x=np.column_stack((a['x'],a['y']));mesh=Voronoi(x);pairs=[];vertices=[]
 for ij,rr in zip(mesh.ridge_points,mesh.ridge_vertices):
  if len(rr)==2 and min(rr)>=0:pairs.append(ij);vertices.append(rr)
 i,j=np.asarray(pairs).T;q=mesh.vertices[np.asarray(vertices)];A=np.linalg.norm(q[:,1]-q[:,0],axis=1);dx=x[j]-x[i];d=np.linalg.norm(dx,axis=1);e=dx/d[:,None];t=q.mean(axis=1)-.5*(x[i]+x[j]);t-=np.sum(t*e,axis=1)[:,None]*e
 def force(P):
  z=-A[:,None]*(.5*(P[i]+P[j])[:,None]*e+((P[j]-P[i])/d)[:,None]*t);F=np.zeros((n,2));np.add.at(F,i,z);np.add.at(F,j,-z);return F
 return force
lp,lm=map(geometry,states);p,m=states;Pp=p['pressure'];Pm=m['pressure'];deltaP=Pp-Pm;meanP=.5*(Pp+Pm)
actual=(lp(Pp)-lm(Pm))/(2*eps)
pressure=.5*(lp(deltaP)+lm(deltaP))/(2*eps)
geometric=(lp(meanP)-lm(meanP))/(2*eps)
iso=[initial['pressure']*(initial['volume']/a['volume'])**gamma for a in states]
isentropic=(lp(iso[0])-lm(iso[1]))/(2*eps)
entropy=(lp(Pp-iso[0])-lm(Pm-iso[1]))/(2*eps)
def departure(a):
 x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']));om=np.sqrt(gravity_coefficient(np.linalg.norm(x,axis=1)))
 return v-om[:,None]*np.column_stack((-x[:,1],x[:,0]))
du=(departure(p)-departure(m))/(2*eps)
dx=(np.column_stack((p['x'],p['y']))-np.column_stack((m['x'],m['y'])))/(2*eps)
r0=np.hypot(ic['x']-3,ic['y']-3);masks={'inner':r0<.5,'shear_band':(r0>=.3)&(r0<.5),'central':r0<.2,'edge':(r0>=.5)&(r0<.6),'core':(r0>=.6)&(r0<1.9)}
sel=r0<1.9;err=float(np.linalg.norm((actual-pressure-geometric)[sel])/np.linalg.norm(actual[sel]));err2=float(np.linalg.norm((actual-isentropic-entropy)[sel])/np.linalg.norm(actual[sel]));assert max(err,err2)<1e-8
parts={'actual':actual,'pressure_value':pressure,'geometry':geometric,'seed_entropy_isentropic':isentropic,'entropy_excess_since_seed':entropy};rows={}
for name,mask in masks.items():
 rows[name]={k:{'departure_response_power':float(np.sum(du[mask]*F[mask])),
 'minus_displacement_dot_force_response':float(-np.sum(dx[mask]*F[mask]))} for k,F in parts.items()}
out={'time':time,'seed_time':t0,'epsilon':eps,'pressure_geometry_split_relative_closure':err,'isentropic_entropy_split_relative_closure':err2,'regions':rows,
 'scope':'Instantaneous reversible-pressure contribution sum delta_u dot delta_F at saved finite +/- trajectories. Symmetric exact algebraic pressure-value/geometry split; second split uses fixed seed-t2 entropy and subsequent excess. Omits TCF/HG/gravity/reference-velocity derivative, so not full norm-growth budget. No eigenmode or negative full-system stiffness claim. Independent finite Voronoi faces interpreted only for interior regions.'}
with (root/f'pressure_response_{args.step:06d}.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))

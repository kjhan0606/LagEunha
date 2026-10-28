"""Gravity-only control from the same coupled seed states; no hydro changes."""
import json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');src=b/'Codex_coupled_onset_seed_20261018'
root=b/'Codex_coupled_orbit_control_20261018';root.mkdir(exist_ok=False)
names=['reference','plus','minus'];initial=[];final=[]
for name in names:
 a,t0=read(src/name/'khout.000184.dat');z,t1=read(src/name/'khout.000224.dat')
 assert np.array_equal(a['id'],z['id']);initial.append(a);final.append(z)
def vectors(states):
 return np.array([np.column_stack((a['x']-3,a['y']-3)) for a in states]),np.array([np.column_stack((a['vx'],a['vy'])) for a in states])
x0,v0=vectors(initial);xh,vh=vectors(final);m=initial[0]['mass'];eps=.001
active=np.linalg.norm(x0,axis=2)<2.5
def rhs(x,v):
 r=np.linalg.norm(x,axis=2);acc=-gravity_coefficient(r)[...,None]*x
 return np.where(active[...,None],v,0),np.where(active[...,None],acc,0)
def evolve(n):
 dt=(t1-t0)/n;x=x0.copy();v=v0.copy()
 for step in range(n):
  a,b1=rhs(x,v);c,d=rhs(x+dt*a/2,v+dt*b1/2);e,f=rhs(x+dt*c/2,v+dt*d/2);g,h=rhs(x+dt*e,v+dt*f)
  x+=dt*(a+2*c+2*e+g)/6;v+=dt*(b1+2*d+2*f+h)/6
 assert np.array_equal(np.linalg.norm(x,axis=2)<2.5,active)
 return x,v
def metrics(x,v):
 om2=gravity_coefficient(np.linalg.norm(x,axis=2));vc=np.sqrt(om2)[...,None]*np.stack((-x[:,:,1],x[:,:,0]),axis=2)
 dx=(x[1]-x[2])/(2*eps);dv=(v[1]-v[2])/(2*eps);du=((v[1]-vc[1])-(v[2]-vc[2]))/(2*eps)
 return {'velocity':float(np.sum(m[:,None]*dv**2)),
  'circular_departure_velocity':float(np.sum(m[:,None]*du**2)),
  'omega_weighted_displacement':float(np.sum(m[:,None]*om2[0,:,None]*dx**2))}
start=metrics(x0,v0);hydro=metrics(xh,vh);cases={};endpoints={}
for n in [200,400]:
 x,v=evolve(n);met=metrics(x,v);cases[str(n)]={'dt':(t1-t0)/n,'metrics':met,'ratios':{k:met[k]/start[k] for k in met}}
 endpoints[n]=(x,v)
x,v=endpoints[400];xc,vc=endpoints[200]
dv=(v[1]-v[2])/(2*eps);dvcoarse=(vc[1]-vc[2])/(2*eps)
response_err=float(np.sqrt(np.sum(m[:,None]*(dv-dvcoarse)**2)/np.sum(m[:,None]*dv**2)))
np.savez(root/'fine_endpoints.npz',x=x,v=v)
out={'start_time':t0,'end_time':t1,'initial_metrics':start,'hydro_ratios':{k:hydro[k]/start[k] for k in hydro},
 'gravity_only':cases,'gravity_timestep_velocity_response_relative_difference':response_err,
 'scope':'Independent RK4 gravity-only evolution of identical initial reference/+/- states with same spline potential; frozen exterior retained. No pressure/AV/HG or force feedback from geometry, no hydrodynamic stability claim. Squared positive metrics include orbital phase/shear effects. Control is not obtained by altering production code.'}
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))

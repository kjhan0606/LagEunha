"""Isentropic pressure proxy along gravity-only control, no pressure feedback."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.spatial import Voronoi
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce');src=b/'Codex_coupled_onset_seed_20261018';root=b/'Codex_coupled_orbit_control_20261018'
end=np.load(root/'fine_endpoints.npz');gamma=float(np.float32(1.666667));eps=.001
initial=[];hydro=[]
for name in ['reference','plus','minus']:
 a,t=read(src/name/'khout.000184.dat');z,t1=read(src/name/'khout.000224.dat');initial.append(a);hydro.append(z)
mask=np.hypot(initial[0]['x']-3,initial[0]['y']-3)<1.9;ids=np.flatnonzero(mask)
def volumes(x):
 mesh=Voronoi(x);v=[]
 for i in ids:
  reg=mesh.regions[mesh.point_region[i]];assert reg and min(reg)>=0
  q=mesh.vertices[reg]-x[i];v.append(abs(np.dot(q[:,0],np.roll(q[:,1],-1))-np.dot(q[:,1],np.roll(q[:,0],-1)))/2)
 return np.array(v)
proxy=[];areas=[];errors=[]
for k,a in enumerate(initial):
 x=np.column_stack((a['x'],a['y']));V0=volumes(x);V1=volumes(end['x'][k]+3)
 error=float(np.max(abs(V0/a['volume'][mask]-1)));assert error<1e-9;errors.append(error)
 proxy.append(a['pressure'][mask]*(V0/V1)**gamma);areas.append(a['volume'][mask]*V1/V0)
def acoustic(p,V):
 dp=(p[1]-p[2])/(2*eps);return float(np.sum(V[0]*dp**2/(gamma*p[0])))
p0=[a['pressure'][mask] for a in initial];V0=[a['volume'][mask] for a in initial]
ph=[a['pressure'][mask] for a in hydro];Vh=[a['volume'][mask] for a in hydro]
n0=acoustic(p0,V0);nh=acoustic(ph,Vh);nc=acoustic(proxy,areas)
out={'region':'Initial reference r<1.9 fixed particle IDs','initial_geometry_relative_errors':errors,
 'initial_acoustic_norm':n0,'hydro_acoustic_ratio':nh/n0,'gravity_geometry_pressure_proxy_ratio':nc/n0,
 'scope':'Pressure proxy P=P_start*(V_start/V_ballistic)^gamma along independent gravity-only trajectories; entropy held fixed in each case. No pressure feedback, no thermodynamic evolution, not a hydrodynamic solution or causal growth fraction.'}
with (root/'pressure_proxy.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))

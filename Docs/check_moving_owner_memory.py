"""Moving periodic 2D auxiliary subsystem only; not Euler or Kepler.

Constant masses and acoustic scale, x_dot=v, m*v_dot=-B(x)^T*z,
z_dot=B(x)*v. Every RK stage rebuilds support and operator.
"""
import json,sys
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_reversible_memory_moving_20261018')
from smooth_owner_memory import build_row
nside=12;n=nside*nside;h=1/nside;m=h*h;c=1.
rng=np.random.default_rng(20261018)
grid=np.array([(i,j) for i in range(nside) for j in range(nside)])
x0=(grid+.5+rng.uniform(-.12,.12,(n,2)))/nside
v0=np.column_stack((.2*np.sin(2*np.pi*x0[:,1]),.02*(-1.)**grid[:,0]))
y0=np.array([x0,v0,np.zeros((n,2))]);initial_energy=.5*m*np.sum(v0*v0)
def evaluate(y,keys=False):
    x,v,z=y;tree=cKDTree(x%1,boxsize=1.)
    force=np.zeros((n,2));dz=np.zeros((n,2));supports=[]
    for i in range(n):
        ids,r,_=build_row(x,i,tree,box=1.)
        b=np.sqrt(m)*c/h*r
        dz[i]=b@v[ids];np.add.at(force,ids,-b[:,None]*z[i])
        if keys:supports.append(tuple(sorted(map(int,ids))))
    return (np.array([v,force/m,dz]),supports) if keys else np.array([v,force/m,dz])
def energy(y):return .5*(m*np.sum(y[1]**2)+np.sum(y[2]**2))
_,initial_support=evaluate(y0,True)
rows=[];endpoints=[]
for dt in [.002,.001,.0005]:
    y=y0.copy()
    for _ in range(round(.02/dt)):
        a=evaluate(y);b=evaluate(y+dt*a/2);cc=evaluate(y+dt*b/2);d=evaluate(y+dt*cc)
        y+=dt/6*(a+2*b+2*cc+d)
    _,final_support=evaluate(y,True)
    rows.append({'dt':dt,'time':.02,'energy_relative_error':float(energy(y)/initial_energy-1),
                 'momentum_absolute_error':float(np.linalg.norm(m*(y[1]-v0).sum(axis=0))),
                 'auxiliary_energy_fraction':float(.5*np.sum(y[2]**2)/initial_energy),
                 'owners_with_changed_support':sum(a!=b for a,b in zip(initial_support,final_support))})
    endpoints.append(y)
    print(json.dumps({'completed':rows[-1]}),file=sys.stderr,flush=True)
differences=[float(np.linalg.norm(endpoints[j]-endpoints[j+1])) for j in [0,1]]
assert all(abs(row['energy_relative_error'])<1e-5 for row in rows)
assert all(row['momentum_absolute_error']<1e-12 for row in rows)
assert differences[1]<differences[0]/8
print(json.dumps({'scope':__doc__,'cells':n,'rows':rows,
                  'successive_full_state_L2_differences':differences,
                  'ratio':differences[0]/differences[1],
                  'limitations':'Constant acoustic coefficient and auxiliary subsystem only; no pressure/gravity/heat, no production restart or parallel implementation.'},indent=2))

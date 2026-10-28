"""Moving periodic auxiliary subsystem with reconstructed spin; not Kepler.

xdot=v, m*vdot=-B^T z, zdot=Bv+W(v,x)z. Nonzero initial z exercises
transport immediately. Rebuild support and velocity gradient every RK stage.
"""
import json,sys
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_memory_circular_work_20261018')
from smooth_owner_memory import build_row

nside=12; n=nside*nside; h=1/nside; m=h*h
rng=np.random.default_rng(20261018)
grid=np.array([(i,j) for i in range(nside) for j in range(nside)])
x0=(grid+.5+rng.uniform(-.12,.12,(n,2)))/nside
v0=np.column_stack((.2*np.sin(2*np.pi*x0[:,1]),.02*(-1.)**grid[:,0]))
z0=.001*rng.normal(size=(n,2))
y0=np.array([x0,v0,z0])
max_power_residual=0.; max_frequency_bound=0.; max_spin=0.

def evaluate(y,keys=False):
    global max_power_residual,max_frequency_bound,max_spin
    x,v,z=y; tree=cKDTree(x%1,boxsize=1.)
    force=np.zeros((n,2)); dz=np.zeros((n,2)); supports=[]
    columns=np.zeros(n); rowmax=0.; spinmax=0.
    for i in range(n):
        ids,res,meta=build_row(x,i,tree,box=1.)
        radius=meta['radius']; dx=x[ids]-x[i]; dx-=np.rint(dx); q=dx/radius
        weights=np.maximum(1-np.sum(q*q,axis=1),0)**4
        P=np.column_stack((np.ones(len(ids)),q,q[:,0]**2,q[:,0]*q[:,1],q[:,1]**2))
        coef=np.linalg.lstsq(np.sqrt(weights)[:,None]*P,
                            np.sqrt(weights)[:,None]*v[ids],rcond=None)[0]
        spin=.5*(coef[1,1]-coef[2,0])/radius
        b=np.sqrt(m)/h*res
        dz[i]=b@v[ids]+spin*np.array([-z[i,1],z[i,0]])
        np.add.at(force,ids,-b[:,None]*z[i])
        entries=abs(b)/np.sqrt(m); rowmax=max(rowmax,float(entries.sum()))
        np.add.at(columns,ids,entries); spinmax=max(spinmax,abs(float(spin)))
        if keys:supports.append(tuple(sorted(map(int,ids))))
    residual=abs(float(np.sum(force*v)+np.sum(z*dz)))
    max_power_residual=max(max_power_residual,residual)
    # Bound for the frozen B,W skew block; not a full nonlinear Jacobian bound.
    max_frequency_bound=max(max_frequency_bound,np.sqrt(rowmax*columns.max())+spinmax)
    max_spin=max(max_spin,spinmax)
    result=np.array([v,force/m,dz])
    return (result,supports) if keys else result

def energy(y):return .5*(m*np.sum(y[1]**2)+np.sum(y[2]**2))
initial_energy=energy(y0)
_,initial_support=evaluate(y0,True)
rows=[]; endpoints=[]
for dt in (.002,.001,.0005):
    y=y0.copy()
    for _ in range(round(.02/dt)):
        a=evaluate(y);b=evaluate(y+dt*a/2);c=evaluate(y+dt*b/2);d=evaluate(y+dt*c)
        y+=dt/6*(a+2*b+2*c+d)
    assert np.isfinite(y).all()
    _,support=evaluate(y,True)
    row={'dt':dt,'time':.02,'energy_relative_error':float(energy(y)/initial_energy-1),
         'momentum_absolute_error':float(np.linalg.norm(m*(y[1]-v0).sum(axis=0))),
         'owners_with_changed_support':sum(a!=b for a,b in zip(initial_support,support))}
    rows.append(row);endpoints.append(y)
    print(json.dumps(row),file=sys.stderr,flush=True)
differences=[float(np.linalg.norm(endpoints[j]-endpoints[j+1])) for j in (0,1)]
assert differences[1]<differences[0]/8
assert max(abs(r['energy_relative_error']) for r in rows)<1e-5
assert max(r['momentum_absolute_error'] for r in rows)<1e-12
assert max_power_residual<1e-12
print(json.dumps({'scope':__doc__,'rows':rows,'successive_full_state_L2_differences':differences,
 'ratio':differences[0]/differences[1],'max_stage_power_residual':max_power_residual,
 'maximum_reconstructed_spin':max_spin,'maximum_frozen_block_frequency_bound':max_frequency_bound,
 'largest_dt_times_frequency_bound':.002*max_frequency_bound,
 'limitations':'No pressure, gravity, heat, variable acoustic coefficient or full Euler. Frozen-block bound omits derivatives of B and W; not a nonlinear stability proof. No production implementation or long-time accuracy claim.'},indent=2))

"""Independent quadrature check of density moments in an evolved checkpoint."""
import importlib.util,json,os,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
sys.path[:0]=[str(base/'Codex_continuous_quadratic_evolution_modules_20261018'),
              str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
import moment_ale_evolve as model
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from run_canonical_av_branches_20261018 import read
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
geom=load('geom',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
clip=load('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py')
root=base/'Codex_continuous_quadratic_evolution_dt1e4_20261018'
saved=np.load(root/'step0001.npz');mass=saved['mass'];ids=saved['ids'];a=saved['state'][:-4].reshape(-1,20)
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
builder=EmbeddedLaguerreGeometry(np.column_stack([ic['x']-3,ic['y']-3]),ids,clip.clip)
polys,_=builder(a[:,:2],a[:,2],float(saved['time']))
os.environ['LAGFORCE_DENSITY_QUADRATURE_ORDER']='12';errors=[];minimum=float('inf');worst=[]
for i,(p,m) in enumerate(zip(polys,mass)):
    q=a[i,3:5];Q=a[i,5:9].reshape(2,2)
    density=model.Density(p,m,q,Q)
    x,dv=geom.quadrature(p,24);rho=density.value(x)
    assert np.isfinite(rho).all() and np.all(rho>0)
    minimum=min(minimum,float(rho.min()));z=(x-density.origin)/density.scale
    F=np.column_stack([np.ones(len(z)),z,z[:,0]**2,z[:,0]*z[:,1],z[:,1]**2])
    reconstructed=(dv*rho)@F/m
    c=q/m;o=density.origin
    centered=Q/m-np.outer(o,c)-np.outer(c,o)+np.outer(o,o)
    target=np.r_[1.,(c-o)/density.scale,centered[0,0]/density.scale**2,
                 centered[0,1]/density.scale**2,centered[1,1]/density.scale**2]
    error=float(abs(reconstructed-target).max());errors.append(error);worst.append((error,int(ids[i])))
report=dict(scope=__doc__,time=float(saved['time']),cells=len(ids),fit_order=12,independent_order=24,
    maximum_dimensionless_moment_error=float(max(errors)),minimum_sampled_density=minimum,
    worst_cells=sorted(worst,reverse=True)[:5],
    limits='Spatial moment fit only, at positive quadrature samples. Does not establish density pointwise accuracy, potential closure, time convergence, or long-term stability.')
print(json.dumps(report,indent=2))
assert max(errors)<1e-9

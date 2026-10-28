"""Read-only admissibility audit of actual-t2 lifted/remapped moments.
No covariance clipping, density repair, or checkpoint replacement.
"""
import hashlib,importlib.util,json,sys,time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
sys.path.insert(0,str(base/'Codex_coupled_relaxation_evolution_20261018'))
sys.path.insert(0,str(base))
import moment_ale_evolve as model
from run_canonical_av_branches_20261018 import read
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
geometry=load('saved_remap_geometry',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
clipper=load('saved_remap_clipper',base/'Codex_weighted_compression_geometry_20261017/audit.py')
snapshot,t=read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat')
x=np.column_stack([snapshot['x']-3,snapshot['y']-3])
weights=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights']
results=[]
for fraction in (.25,1.):
    start=time.monotonic();path=base/f'Codex_checkpoint_weight_remap_20261018/state_w{fraction:g}_q8.npz'
    data=np.load(path);ids=data['ids'];q=data['state'];polys=geometry.polygons(x,fraction*weights,ids,clipper.clip)
    errors=[];failures=[];negative=[];min_eigen=float('inf');min_thermal=float('inf');max_fit=0.;min_density=float('inf')
    for local,(index,p,a) in enumerate(zip(ids,polys,q)):
        # Stored 17 uses symmetric Q/Pi and includes potential energy; research raw18 uses full matrices.
        raw=np.r_[a[0:3],a[3],a[4],a[4],a[5],a[6:12],a[12],a[13],a[13],a[14],a[15]]
        assert len(raw)==18
        errors.append(abs(geometry.area(p)/data['volume'][local]-1))
        try:
            state=model.decode(p,raw)
            eig=np.linalg.eigvalsh(state['R']);min_eigen=min(min_eigen,float(eig.min()))
            if eig.min()<0:negative.append(dict(id=int(index),minimum=float(eig.min()),maximum=float(eig.max())))
            min_thermal=min(min_thermal,float(a[15]-.5*(a[12]+a[14])))
            max_fit=max(max_fit,state['density'].residual)
            min_density=min(min_density,float(state['density'].value(p).min()))
        except Exception as e:failures.append(dict(id=int(index),error=repr(e)))
    results.append(dict(fraction=fraction,source_time=t,source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        cells=len(ids),geometry_max_relative_area_error=max(errors),density_decode_failures=failures,
        negative_residual_covariance=negative,minimum_residual_eigenvalue=min_eigen,
        minimum_thermal_energy=min_thermal,maximum_density_fit_residual=max_fit,minimum_vertex_density=min_density,seconds=time.monotonic()-start))
    assert max(errors)<1e-7
print(json.dumps(dict(scope=__doc__,results=results,
    limits='Successful density moment fit and sampled vertex density do not establish complete joint spatial-velocity realizability. Negative raw covariance may be subtraction roundoff; this audit does not repair it.'),indent=2))

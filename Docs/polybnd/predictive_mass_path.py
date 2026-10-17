"""Mass-only evaluation of an existing dense path, without re-solving weights."""
import numpy as np
import moment_ale_evolve as model
from kinetic_factor_rhs import ALPHA
from affine_kinetic_mass import make_bins,make_bin_arrays,arrays,evaluate


def evaluate_path(t,y,L,mass,geometry_rates,provider):
    a=y[:-4].reshape(-1,20);polys,edges=model.geometry(a[:,:2],a[:,2],t)
    raw=np.column_stack([mass,a[:,3:]]);states=[model.decode(v,r) for v,r in zip(polys,raw)]
    material=np.array([model.velocity(s,x) for s,x in zip(states,a[:,:2])]);gen,z=geometry_rates
    trace,closure=provider(polys,edges,raw)
    if provider.vector_geometry:
        data=make_bin_arrays(edges,a[:,:2],states,L,mass,gen,trace,provider.segments,model.boundary_velocity,ALPHA)
    else:
        bins=make_bins(edges,a[:,:2],states,L,mass,gen,trace,provider.segments,model.boundary_velocity,model.velocity,ALPHA)
        data=arrays(bins)
    residual,_=evaluate(data,len(mass),z)
    return -residual,dict(min_pressure=min(s['P'] for s in states),
        generator_defect=float(np.max(abs(gen-material))),
        min_factor_R=min(float(np.linalg.eigvalsh(f@f.T).min()/m) for f,m in zip(L,mass)),**closure)

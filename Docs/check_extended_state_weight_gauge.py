"""Common-weight gauge invariance with retained nonzero covariance state.
Not a proof of persistent equal-weight/rate fixed-mass Voronoi reduction.
"""
import json,os,sys
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_coupled_relaxation_evolution_20261018')
from adaptive_moment_batch import setup
from factor_moment_evolve import rebuild
from kinetic_factor_rhs import rhs
os.environ.update(LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',
    LAGFORCE_VELOCITY_STENCIL='area',LAGFORCE_QUADRATIC_PRESSURE='0',LAGFORCE_KINETIC_SMOOTHING='0',LAGFORCE_VECTOR_KINETIC_ROOT='0',LAGFORCE_PRESSURE_RIEMANN='hllc_davis')
mass,y=setup('shear',4,0.)
factors=[np.sqrt(m)*np.array([[.01,0.],[.003,.007]]) for m in mass]
# Preserve thermal energy when introducing the diagnostic residual energy.
a=y[:-4].reshape(-1,20);a[:,19]+=[.5*np.sum(L*L) for L in factors]
y=rebuild(y,mass,factors);source=y.copy();reference=None;rows=[]
for offset in (0.,.125,-.125):
    state=source.copy();state[:-4].reshape(-1,20)[:,2]+=offset
    d,diag=rhs(0.,state,mass,covariance_factors=factors,return_covariance_terms=True)
    if reference is None:reference=d.copy()
    error=float(np.max(abs(d-reference)));assert error<1e-10
    assert np.array_equal(state[:-4].reshape(-1,20)[:,3:],source[:-4].reshape(-1,20)[:,3:])
    rates=d[:-4].reshape(-1,20)
    rows.append(dict(offset=offset,maximum_rhs_difference=error,mass_flux=diag['mass_flux'],
        relative_weight_rate_range=float(np.ptp(rates[:,2])),generator_velocity_defect=diag['generator_velocity_defect']))
print(json.dumps(dict(scope=__doc__,rows=rows,residual_energy=float(sum(.5*np.sum(L*L) for L in factors)),
    conclusion='No moment or covariance reset. Gauge invariance is necessary but weaker than the approved persistent equal-weight/equal-rate reduction.'),indent=2))

"""Check gravity-to-factor tangent against saved actual quadratic cell sources."""
import json
from pathlib import Path
import numpy as np
from moment_gravity_source import covariance_source,covariance_generator

base=Path('/gpfs/kjhan/LagForce')
h=np.load(base/'Codex_continuous_initial_rhs_quadratic_split_20261018/rhs_r0.55_f8.npz')
g=np.load(base/'Codex_continuous_quadratic_gravity_20261018/gravity_q24.npz')
y=h['state'][:-4].reshape(-1,20);raw=np.column_stack([h['mass'],y[:,3:]])
errors=[];rates=[];generators=[]
for r,d,L in zip(raw,g['source'],h['covariance_factors']):
    G=covariance_source(r,d);A=covariance_generator(L,G);R=L@L.T
    errors.append(np.max(abs(A@R+R@A.T-G)));rates.append(G);generators.append(A)
rates=np.asarray(rates);generators=np.asarray(generators)
error=float(max(errors));scale=float(abs(rates).max())
assert error<1e-12*scale
rank_guard=False
try:covariance_generator(np.zeros((2,2)),np.eye(2))
except ValueError:rank_guard=True
assert rank_guard
report=dict(scope=__doc__,cells=len(raw),maximum_covariance_source=scale,
    maximum_tangent_identity_error=error,maximum_generator=float(abs(generators).max()),
    rank_deficient_guard=rank_guard,
    limits='Strictly positive covariance only; tangent identity, not finite-step energy, order, moving mesh or long-term stability proof.')
np.savez('gravity_generators.npz',ids=h['ids'],rates=rates,generators=generators)
Path('report.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))

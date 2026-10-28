"""Quadratic cuts versus known roots and analytic positive-part integrals."""
import json,sys
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_quadratic_split_modules_20261018')
from quadratic_upwind_face_nodes import polynomial_cuts
rng=np.random.default_rng(7013);cases=[]
for _ in range(500):
    roots=rng.uniform(-.5,1.5,2);scale=rng.choice([-1.,1.])*10**rng.uniform(-10,10)
    cases.append((scale*np.array([roots.prod(),-roots.sum(),1.]),roots))
cases.extend([(np.array([.25,-1.,1.]),[.5,.5]),(np.array([0.,-1.,1.]),[0.,1.]),(np.array([-.5,1.,0.]),[.5]),(np.zeros(3),[])])
gx,gw=np.polynomial.legendre.leggauss(4);gx=(gx+1)/2;gw/=2;error=0.;root_error=0.
for p,known in cases:
    f=lambda t:p[0]+p[1]*t+p[2]*t*t
    cuts=polynomial_cuts(np.array([f(0)]),np.array([f(.5)]),np.array([f(1)]))
    refcuts=np.unique([0.,1.]+[float(r) for r in known if 0<r<1])
    primitive=lambda t:p[0]*t+p[1]*t*t/2+p[2]*t**3/3
    ref=sum(primitive(b)-primitive(a) for a,b in zip(refcuts[:-1],refcuts[1:]) if f((a+b)/2)>0)
    val=sum((b-a)*np.dot(gw,np.maximum(f(a+(b-a)*gx),0)) for a,b in zip(cuts[:-1],cuts[1:]))
    error=max(error,abs(val-ref)/max(1.,np.linalg.norm(p)))
    for r in refcuts[1:-1]:root_error=max(root_error,float(np.min(abs(cuts-r))))
assert error<1e-13 and root_error<1e-10
print(json.dumps(dict(cases=len(cases),maximum_scaled_integral_error=error,maximum_root_error=root_error,
    scope='Known factorized roots incl two crossings, tangency, endpoints, linear and zero cases. No hydro evolution.'),indent=2))

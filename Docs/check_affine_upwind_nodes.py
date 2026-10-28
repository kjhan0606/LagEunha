"""Independent positive-part integral check for affine upwind cuts."""
import json,sys
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_coupled_relaxation_evolution_20261018')
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_split_face_modules_20261018')
from affine_upwind_face_nodes import AffineUpwindNodes
from kinetic_factor_rhs import ALPHA
rng=np.random.default_rng(8091);maximum=0.;minimum_weight=1.;split=0
def integral_positive(a,b):
    if a>=0 and b>=0:return .5*(a+b)
    if a<=0 and b<=0:return 0.
    if a>0:return .5*a*a/(a-b)
    return .5*b*b/(b-a)
for k in range(300):
    angle=rng.uniform(-np.pi,np.pi);n=np.array([np.cos(angle),np.sin(angle)]);t=np.array([-n[1],n[0]])
    x=np.array([-.5*n,.5*n]);a=-.5*t;b=.5*t
    states=[dict(v=rng.normal(size=2),B=rng.normal(size=(2,2)),c=p.copy()) for p in x]
    mass=np.exp(rng.normal(size=2));L=[rng.normal(size=(2,2))*.1 for _ in range(2)]
    z=rng.normal(size=2);gen=rng.normal(size=(2,2))
    provider=AffineUpwindNodes(x,mass,L,z,4,lambda p:np.zeros(2));nodes=provider(0,(0,1,a,b,n,1.),states,gen)
    def speeds(point,side):
        s=states[side];vel=s['v']+s['B']@(point-s['c'])+ALPHA@L[side].T/np.sqrt(mass[side])
        wn=(n@gen[1]-(gen[1]-gen[0])@(point-x[0]))+(z[0]-z[1])/2
        return vel@n-wn
    for side in (0,1):
        reference=sum(integral_positive(u,v) for u,v in zip(speeds(a,side),speeds(b,side)))/4
        actual=sum(w*np.maximum(speeds(a+g*(b-a),side),0).sum()/4 for g,w in nodes)
        maximum=max(maximum,abs(actual-reference))
    assert abs(sum(w for g,w in nodes)-1)<1e-14
    minimum_weight=min(minimum_weight,min(w for g,w in nodes));split+=provider.split_faces
assert maximum<1e-13 and minimum_weight>0
print(json.dumps(dict(cases=300,integrals=600,maximum_absolute_error=maximum,minimum_weight=minimum_weight,split_faces=split,
    scope='Constant-density affine positive-part integral against closed-form reference. No hydro evolution or smooth-density convergence claim.'),indent=2))

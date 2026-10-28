"""Explicit cut/tangent/flat/quadratic/mirror scalar-versus-batch fixtures.

Synthetic algebra fixtures only: no density fit, RHS or provider evaluation.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--modules',type=Path,required=True)
args=parser.parse_args()
sys.path[:0]=[str(args.modules),str(base/'Codex_compiled_midpoint_modules_20261018'),
    str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),
    str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from contact_prior_density import ContactPriorDensity
from quadratic_velocity_reconstruction import velocity_samples_and_mean,velocity_samples_and_mean_batch
from kinetic_factor_rhs import ALPHA
rng=np.random.default_rng(20261018)
eps=np.finfo(float).eps
results=[]


def comparison(actual,reference):
    assert np.isfinite(actual).all() and np.isfinite(reference).all()
    error=abs(actual-reference)
    bound=64*eps*(1+abs(reference))
    row=dict(bitwise_equal=actual.tobytes()==reference.tobytes(),maximum_abs=float(error.max(initial=0.)),
             maximum_roundoff_bound_fraction=float((error/bound).max(initial=0.)))
    assert np.all(error<=bound),row
    return row


def density(meta):
    d=ContactPriorDensity.__new__(ContactPriorDensity)
    d.mass=1.7;d.origin=np.array([.13,-.2]);d.scale=.71
    d.theta=np.array([.02,-.03,.01,.015,-.01]);d.logZ=.47;d.cut_metadata=meta
    return d


points=np.vstack([rng.uniform(-.6,.6,(127,2)),[[.25,0.],[np.nextafter(.25,-np.inf),0.],[np.nextafter(.25,np.inf),0.]]])
directions=rng.normal(size=points.shape)
directions[-3:]=[[1.,0.],[-1.,0.],[1.,0.]]
for meta in (dict(flat=True,value=3.),dict(flat=False,normal=[1.,0.],offset=.25,low=.01,high=1.),
             dict(flat=False,normal=[.6,.8],offset=.125,low=.01,high=1.)):
    d=density(meta)
    scalar=np.array([d.value_trace(p,n)[0] for p,n in zip(points,directions)])
    actual,branches=d.value_trace_batch(points,directions,return_branch=True)
    expected=[]
    for p,n in zip(points,directions):
        if meta['flat']:expected.append(0);continue
        normal=np.asarray(meta['normal']);signed=float((np.atleast_2d(p)@normal)[0]-meta['offset'])
        derivative=float((np.atleast_2d(n)@normal)[0])
        expected.append(1 if (signed>0 or (signed==0 and derivative>0)) else -1)
    assert np.array_equal(branches,np.array(expected,dtype=np.int8))
    results.append(dict(kind='density',flat=meta['flat'],branches_exact=True,values=comparison(actual,scalar)))
d=density(dict(flat=False,normal=[1.,0.],offset=.25,low=.01,high=1.))
errors=[]
for fn in (d.value_trace,d.value_trace_batch):
    try:fn(np.array([[.25,0.]]),np.array([[0.,1.]]))
    except ValueError as error:errors.append(str(error))
assert len(errors)==2 and errors[0]==errors[1]
results.append(dict(kind='on_cut_tangent',both_raise_same_error=errors[0]))
d=density(dict(flat=True,value=3.))
results.append(dict(kind='flat_ignores_direction_as_scalar',values=comparison(
    d.value_trace_batch(points,np.full_like(points,np.nan)),np.array([d.value_trace(p,[np.nan,np.nan])[0] for p in points]))))

L=np.array([[.02,0.],[.003,.01]])
state=dict(v=np.array([.4,-.7]),c=np.array([.1,-.12]),B=np.array([[.2,-.9],[.6,.13]]))
for quadratic in (False,True):
    s=state.copy()
    if quadratic:
        s['velocity_quadratic']=dict(center=np.array([.12,-.13]),scale=.4,
            projection=rng.normal(size=(3,3))*.03,standard_coefficient=rng.normal(size=(2,3))*.02,
            residual_map=np.array([[.8,.1],[-.1,.9]]),coefficient=rng.normal(size=(2,3))*.01)
    scalar=[velocity_samples_and_mean(s,L,1.7,p,ALPHA) for p in points]
    reference=(np.concatenate([v[0] for v in scalar]),np.concatenate([v[1] for v in scalar]),np.array([v[2] for v in scalar]))
    actual=velocity_samples_and_mean_batch(s,L,1.7,points,ALPHA)
    checks={name:comparison(v,r) for name,v,r in zip(('samples','xi','shared_mean'),actual,reference)}
    normals=directions/np.linalg.norm(directions,axis=1)[:,None]
    wn=np.linspace(-.1,.1,len(points))
    mirror=np.array([v-2*((v@n-w)[:,None]*n) for v,n,w in zip(actual[0],normals,wn)])
    old=np.array([v-2*((v@n-w)[:,None]*n) for v,n,w in zip(reference[0],normals,wn)])
    checks['mirror_samples']=comparison(mirror,old)
    results.append(dict(kind='velocity',quadratic=quadratic,checks=checks))
report=dict(scope=__doc__,status='passed',fixtures=results,
    source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                   (Path(__file__),args.modules/'contact_prior_density.py',args.modules/'quadratic_velocity_reconstruction.py',args.modules/'kinetic_factor_rhs.py')},
    limits='Algebra fixtures, not saved-stage/provider/long-time validation. Exact cut branches; numerical bound64 eps*(1+abs(scalar)) is a local fixture gate only.')
print(json.dumps(report,indent=2))

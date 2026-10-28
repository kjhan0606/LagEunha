"""Same saved stage0/half nodes: scalar/batch fields and exact cut branches.

Reuses archived fitted theta and old polygons; no fit, RHS or provider call.
"""
import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--modules',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
sys.path[:0]=[str(args.modules),str(base/'Codex_compiled_midpoint_modules_20261018'),
    str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),
    str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from saved_contact_density_restore import restore_density
import moment_ale_evolve as model
import quadratic_velocity_reconstruction as quad
from kinetic_factor_rhs import ALPHA
os.environ.update(LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',
    LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area')
model.GAMMA=float(np.float32(1.666667))
args.output.mkdir(exist_ok=False)
(args.output/'manifest.json').write_text(json.dumps(dict(scope=__doc__,modules=str(args.modules),
    source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        (Path(__file__),args.modules/'contact_prior_density.py',args.modules/'quadratic_velocity_reconstruction.py',
         args.modules/'saved_contact_density_restore.py')},provider_evaluations=0,newton_fits=0),indent=2))
eps=np.finfo(float).eps
reports=[]


def compare(actual,old):
    assert actual.shape==old.shape and np.isfinite(actual).all() and np.isfinite(old).all()
    delta=abs(actual-old);bound=64*eps*(1+abs(old))
    result=dict(bitwise_equal=actual.tobytes()==old.tobytes(),maximum_abs=float(delta.max(initial=0)),
        maximum_roundoff_bound_fraction=float((delta/bound).max(initial=0)))
    assert np.all(delta<=bound),result
    return result


for label in ('stage0','midpoint'):
    started=time.perf_counter()
    root=base/'Codex_adjacency_prior_impact_20261018'/label
    common=np.load(root/'common_inputs.npz')
    priors=json.loads((root/'priors.json').read_text())['old']
    fits=json.loads((root/'fit_residuals.json').read_text())
    saved=np.load(base/'Codex_compiled_midpoint_once_20261018'/label/'inputs.npz')
    assert np.array_equal(common['state'],saved['state']) and np.array_equal(common['ids'],saved['ids'])
    raw=common['raw'];mass=saved['mass'];L=saved['factors'];ids=saved['ids']
    offsets=common['polygon_offsets'];vertices=common['polygon_vertices']
    polygons=[vertices[a:b] for a,b in zip(offsets[:-1],offsets[1:])]
    densities=[restore_density(p,r,meta,old['old']) for p,r,meta,old in zip(polygons,raw,priors,fits)]
    gravity_offsets=saved['gravity_offsets'];gravity_nodes=saved['gravity_nodes']
    for i,d in enumerate(densities):
        lo,hi=gravity_offsets[i:i+2];old=gravity_nodes[lo:hi]
        assert np.array_equal(d.points,old[:,:2]) and np.array_equal(d.weighted_mass,old[:,2])
    states=[model.decode(p,r,density=d) for p,r,d in zip(polygons,raw,densities)]
    for state,factor in zip(states,L):state['R']=factor@factor.T
    # Actual radial_fixed_mean prepare does not use its edge argument.
    quad.prepare(states,polygons,[],mass,L)
    nodes=saved['node_inputs'];owners=saved['owners'];extra=saved['extra']
    count=len(nodes);entries=[[] for _ in ids]
    for k,(i,j) in enumerate(owners):
        entries[i].append((k,0))
        if j>=0:entries[j].append((k,1))
    summaries={};branches_exact=True
    maximum={name:dict(maximum_abs=0.,maximum_roundoff_bound_fraction=0.,bitwise_equal=True)
             for name in ('density','samples','xi','shared_mean','saved_samples','saved_xi')}
    scalar_calls=0;batch_calls=0;sizes=[]
    for i,items in enumerate(entries):
        if not items:continue
        index=np.asarray(items);rows=index[:,0];sides=index[:,1]
        points=nodes[rows,:2];inward=(2*sides-1)[:,None]*nodes[rows,2:4]
        scalar_rho=[];scalar_v=[];scalar_xi=[];scalar_mean=[];scalar_branch=[]
        d=densities[i];meta=d.cut_metadata
        for p,n in zip(points,inward):
            scalar_rho.append(d.value_trace(p,n)[0])
            v,xi,mean=quad.velocity_samples_and_mean(states[i],L[i],mass[i],p,ALPHA)
            scalar_v.append(v[0]);scalar_xi.append(xi[0]);scalar_mean.append(mean)
            if meta['flat']:scalar_branch.append(0)
            else:
                normal=np.asarray(meta['normal'])
                signed=(np.atleast_2d(p)@normal)[0]-meta['offset'];derivative=(np.atleast_2d(n)@normal)[0]
                assert signed!=0 or derivative!=0
                scalar_branch.append(1 if (signed>0 or (signed==0 and derivative>0)) else -1)
        rho,branch=d.value_trace_batch(points,inward,return_branch=True)
        v,xi,mean=quad.velocity_samples_and_mean_batch(states[i],L[i],mass[i],points,ALPHA)
        assert np.array_equal(branch,np.asarray(scalar_branch,dtype=np.int8))
        old_v=np.array([nodes[k,5+8*side:13+8*side].reshape(4,2) for k,side in items])
        old_xi=np.array([extra[k,8*side:8+8*side].reshape(4,2) for k,side in items])
        checks=dict(density=compare(rho,np.array(scalar_rho)),samples=compare(v,np.array(scalar_v)),
            xi=compare(xi,np.array(scalar_xi)),shared_mean=compare(mean,np.array(scalar_mean)),
            saved_samples=compare(v,old_v),saved_xi=compare(xi,old_xi))
        for name,row in checks.items():
            maximum[name]['maximum_abs']=max(maximum[name]['maximum_abs'],row['maximum_abs'])
            maximum[name]['maximum_roundoff_bound_fraction']=max(maximum[name]['maximum_roundoff_bound_fraction'],row['maximum_roundoff_bound_fraction'])
            maximum[name]['bitwise_equal'] &= row['bitwise_equal']
        scalar_calls+=len(points);batch_calls+=1;sizes.append(len(points))
        if (i+1)%200==0:print(json.dumps(dict(stage=label,cells=i+1)),flush=True)
    report=dict(stage=label,time=float(saved['time']),nodes=count,cells=len(ids),scalar_side_calls=scalar_calls,
        cell_batches=batch_calls,maximum_cell_batch=max(sizes),branches_exact=True,density_newton_fits=0,
        stored_volume_points_and_mass_weights_bitwise=True,comparisons=maximum,seconds=time.perf_counter()-started,
        limits='Scalar/batch fields at stored converged nodes only. No massroot, RHS or provider. Face splitting/order and complete integration tested separately.')
    (args.output/(label+'.json')).write_text(json.dumps(report,indent=2));reports.append(report)
out=dict(status='passed',scope=__doc__,stages=reports,
    source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                   (Path(__file__),args.modules/'contact_prior_density.py',args.modules/'quadratic_velocity_reconstruction.py',args.modules/'saved_contact_density_restore.py')})
(args.output/'report.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)

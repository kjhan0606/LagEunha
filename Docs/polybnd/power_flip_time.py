"""Candidate copower-circle times on polynomial site/weight paths.

Four copower-circular sites need not form a feasible mesh event. Callers must
verify adjacency and actual geometry; this module never changes a trajectory.
"""
import itertools
import numpy as np
from numpy.polynomial import Polynomial


def determinant(paths):
    paths = np.asarray(paths,dtype=float)
    if paths.ndim!=3 or paths.shape[:2]!=(4,3) or not np.all(np.isfinite(paths)):
        raise ValueError('four finite x,y,weight polynomial paths required')
    rows = []
    for site in paths[1:]:
        dx,dy,dw = (Polynomial(v) for v in site-paths[0])
        rows.append((dx,dy,dx*dx+dy*dy-dw))
    result = Polynomial([0.])
    for permutation in itertools.permutations(range(3)):
        inversions = sum(permutation[i]>permutation[j] for i in range(3) for j in range(i+1,3))
        term = Polynomial([(-1.)**inversions])
        for row,col in enumerate(permutation):
            term *= rows[row][col]
        result += term
    return result


def candidate_times(paths):
    polynomial = determinant(paths)
    scale = np.max(abs(polynomial.coef))
    if scale==0:
        raise ValueError('persistent copower degeneracy requires separate handling')
    roots = (polynomial/scale).roots()
    return sorted(float(z.real) for z in roots
                  if abs(z.imag)<1e-9*(1+abs(z.real)) and 0<z.real<=1)


def disappearing_edge_candidates(paths,start_edges,end_edges):
    """Propose interior flips visible in an endpoint adjacency comparison.

    This deliberately reports unresolved edges. Boundary events, transient
    flips returning to the original adjacency, and persistent degeneracies
    require other handling. Returned roots still need actual geometry checks.
    """
    paths=np.asarray(paths,dtype=float)
    if paths.ndim!=3 or paths.shape[1]!=3 or not np.all(np.isfinite(paths)):
        raise ValueError('finite site x,y,weight polynomial paths required')
    def pairs(edges):
        result=set()
        for edge in edges:
            i,j=map(int,edge[:2])
            if i<0 or j<0:continue
            if i==j or max(i,j)>=len(paths):raise ValueError('invalid interior edge')
            result.add(tuple(sorted((i,j))))
        return result
    start,end=pairs(start_edges),pairs(end_edges)
    neighbors=[set() for _ in paths]
    for i,j in start:neighbors[i].add(j);neighbors[j].add(i)
    candidates=[];unresolved=[]
    for i,j in sorted(start-end):
        common=sorted(neighbors[i]&neighbors[j])
        if len(common)!=2:
            unresolved.append(dict(edge=[i,j],reason='ambiguous adjacent sites',common=common))
            continue
        ids=[i,j,*common]
        try:roots=candidate_times(paths[ids])
        except ValueError as exc:
            unresolved.append(dict(edge=[i,j],reason=str(exc)));continue
        if not roots:unresolved.append(dict(edge=[i,j],reason='no real event root'))
        candidates.extend(dict(edge=[i,j],sites=ids,phase=s) for s in roots)
    return dict(candidates=sorted(candidates,key=lambda row:row['phase']),unresolved=unresolved)


def rejected_step_proposal(record,edge_geometry,minimum_dt=0.):
    """Suggest a shorter retry; every ordinary acceptance gate still applies.

    edge_geometry(x,w,t) must return the actual edges including their lengths.
    Only the unmodified quadratic midpoint site/weight polynomial is supported.
    """
    if 'geometry_lift' in record:raise ValueError('event proposal requires quadratic geometry')
    dt=float(record['dt'])
    if not np.isfinite(dt) or dt<=0 or not np.isfinite(minimum_dt) or minimum_dt<0:
        raise ValueError('positive finite step and nonnegative minimum required')
    y,d0,dm=(np.asarray(record[k],dtype=float)[:-4].reshape(-1,20) for k in ('y','d0','dm'))
    paths=np.stack([y[:,:3],dt*d0[:,:3],dt*(dm-d0)[:,:3]],axis=-1)
    def edges(s):
        a=paths@np.array([1.,s,s*s])
        return edge_geometry(a[:,:2],a[:,2],record['t']+s*dt)
    result=disappearing_edge_candidates(paths,edges(0.),edges(1.))
    result.update(dt=dt/2,selection='halving')
    for candidate in result['candidates']:
        s=candidate['phase']
        # Require meaningful contraction, and avoid proposing an exhausted step.
        if not minimum_dt<=dt*s or not 0<s<=.9:continue
        delta=min(.001,s/2,(1-s)/2)
        lengths=[sum(float(e[-1]) for e in edges(phase)
                     if set(e[:2])==set(candidate['edge'])) for phase in (s-delta,s+delta)]
        candidate['face_lengths_before_after']=lengths
        if lengths[0]>0 and lengths[1]==0:
            result.update(dt=dt*s,selection='verified_interior_event',selected=candidate)
            break
    return result

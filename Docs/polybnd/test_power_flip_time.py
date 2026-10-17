import numpy as np
import pytest
from power_flip_time import determinant,candidate_times
from power_flip_time import disappearing_edge_candidates
from power_flip_time import rejected_step_proposal


def test_known_square_event():
    p=np.zeros((4,3,3));p[:,:2,0]=[[0,0],[1,0],[0,1],[1,.6]];p[3,1,1]=1.
    np.testing.assert_allclose(candidate_times(p),[.4],atol=1e-14)
    q=p.copy();q[:,0]+=np.array([2.,.3,-.2]);q[:,1]+=np.array([-3.,.1,.4]);q[:,2]+=np.array([.7,-.2,.1])
    np.testing.assert_allclose(candidate_times(q),[.4],atol=1e-13)


def test_polynomial_matches_centered_direct_determinant():
    rng=np.random.default_rng(214);p=rng.normal(size=(4,3,3));f=determinant(p)
    for t in np.linspace(0,1,11):
        a=p@np.array([1.,t,t*t]);delta=a[1:]-a[0]
        direct=np.linalg.det(np.column_stack([delta[:,:2],np.sum(delta[:,:2]**2,axis=1)-delta[:,2]]))
        np.testing.assert_allclose(f(t),direct,rtol=1e-11,atol=1e-11)


def test_persistent_degeneracy_is_not_an_isolated_event():
    p=np.zeros((4,3,1));p[:,:2,0]=[[0,0],[1,0],[0,1],[1,1]]
    with pytest.raises(ValueError,match='persistent'):candidate_times(p)


def test_adjacency_selects_disappearing_diagonal():
    p=np.zeros((4,3,3));p[:,:2,0]=[[0,0],[1,0],[0,1],[1,.6]];p[3,1,1]=1.
    boundary=[(0,1),(1,3),(3,2),(2,0)]
    result=disappearing_edge_candidates(p,boundary+[(3,0)],boundary+[(1,2)])
    assert not result['unresolved']
    assert len(result['candidates'])==1
    assert result['candidates'][0]['edge']==[0,3]
    np.testing.assert_allclose(result['candidates'][0]['phase'],.4,atol=1e-14)
    assert disappearing_edge_candidates(p,boundary,boundary)['candidates']==[]


def test_ambiguous_or_degenerate_adjacency_is_reported():
    p=np.zeros((4,3,1));p[:,:2,0]=[[0,0],[1,0],[0,1],[1,1]]
    result=disappearing_edge_candidates(p,[(0,3)],[])
    assert result['unresolved'][0]['reason']=='ambiguous adjacent sites'
    edges=[(0,1),(1,3),(3,2),(2,0)]
    result=disappearing_edge_candidates(p,edges+[(0,3)],edges+[(1,2)])
    assert 'persistent' in result['unresolved'][0]['reason']


def test_retry_proposal_keeps_halving_fallback_and_checks_geometry():
    y=np.zeros((4,20));y[:,:2]=[[0,0],[1,0],[0,1],[1,.6]]
    d=np.zeros_like(y);d[3,1]=1.
    record=dict(t=0.,dt=1.,y=np.r_[y.ravel(),np.zeros(4)],
                d0=np.r_[d.ravel(),np.zeros(4)],dm=np.r_[d.ravel(),np.zeros(4)])
    def edges(x,w,t):
        boundary=[(0,1,1.),(1,3,1.),(3,2,1.),(2,0,1.)]
        return boundary+([(0,3,.4-t)] if t<.4 else [(1,2,t-.4)])
    result=rejected_step_proposal(record,edges)
    assert result['selection']=='verified_interior_event'
    np.testing.assert_allclose(result['dt'],.4)
    assert rejected_step_proposal(record,edges,minimum_dt=.45)['dt']==.5
    assert rejected_step_proposal(record,lambda x,w,t:edges(x,w,0.))['selection']=='halving'
    with pytest.raises(ValueError,match='quadratic'):
        rejected_step_proposal(dict(record,geometry_lift={}),edges)

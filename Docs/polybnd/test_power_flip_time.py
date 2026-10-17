import numpy as np
import pytest
from power_flip_time import determinant,candidate_times


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

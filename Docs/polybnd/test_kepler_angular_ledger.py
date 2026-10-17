import numpy as np
import pytest
from kepler_angular_ledger import components,balance


def fixture():
    a=np.zeros((2,20));a[:,3:5]=[2.,3.];a[:,9:11]=[4.,5.]
    a[:,12]=7.;a[:,13]=3.
    return np.r_[a.ravel(),np.zeros(4)],np.array([2.,2.])


def test_intrinsic_moment_and_active_selection():
    y,m=fixture();c=components(y,m,[1.,1.],[True,False])
    assert c==dict(total=3.,orbital_centroid=-2.,intrinsic=5.)
    assert components(y,m,[1.,1.],[True,True])['total']==6.


def test_boundary_torque_is_shifted_with_momentum():
    y,m=fixture();z=y.copy();a=z[:-4].reshape(-1,20)
    a[0,9:11]+=[.2,.3];a[0,12]+=.7;z[-4:]=[.2,.3,0.,.7]
    for center in ([0.,0.],[1.,1.],[3.,3.]):
        assert abs(balance(y,z,m,center,[True,False])['ledger_residual'])<1e-14


def test_zero_initial_angular_does_not_invent_relative_scale():
    y,m=fixture();y[:40]=0.
    assert balance(y,y,m,[0.,0.],[True,True])['relative_ledger_residual'] is None


def test_invalid_mass_and_shape_rejected():
    y,m=fixture()
    with pytest.raises(ValueError):components(y,[0.,2.],[0.,0.],[True,True])
    with pytest.raises(ValueError):components(y,m,[0.],[True,True])

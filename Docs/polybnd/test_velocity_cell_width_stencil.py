import numpy as np
import pytest
from quadratic_velocity_reconstruction import velocity_tensor_stencil,curvature_fit,stencil_mode
from physical_field_accuracy_audit import environment
from pressure_fit_stencil import tensor_stencil


def test_default_stencil_is_identical(monkeypatch):
    monkeypatch.delenv('LAGFORCE_VELOCITY_STENCIL',raising=False)
    centers=np.array([[0.,0.],[.5,1.],[1.,2.]])
    actual=velocity_tensor_stencil(centers,np.zeros(2),.5,None)
    expected=tensor_stencil(centers,np.zeros(2),.5,3.)
    for a,b in zip(actual,expected):np.testing.assert_array_equal(a,b)


def test_compressed_wall_resolves_planar_curvature(monkeypatch):
    monkeypatch.setenv('LAGFORCE_VELOCITY_STENCIL','cell_width')
    widths=np.array([.25,1.]);x,y=np.meshgrid(np.arange(-4,5)*widths[0],np.arange(4)*widths[1])
    centers=np.column_stack([x.ravel(),y.ravel()]);poly=np.array([[-.125,-.5],[.125,-.5],[.125,.5],[-.125,.5]])
    ids,w=velocity_tensor_stencil(centers,np.zeros(2),.5,poly)
    z=centers[ids]/.5;A=np.column_stack([z[:,0]**2-.1,z[:,0]*z[:,1],z[:,1]**2-.1])
    target=np.column_stack([.2+.3*z[:,0]+z[:,0]**2,np.zeros(len(z))])
    fit=curvature_fit(A,target,w,'tensor_free_affine',zmean=z)
    np.testing.assert_allclose(fit[:,0],[1.,0.,0.],atol=1e-8)
    matrix=np.sqrt(w)[:,None]*np.column_stack([np.ones(len(z)),z,A])
    assert np.linalg.cond(np.vstack([matrix,1e-5*np.eye(6)]))<100
    assert len(np.unique(centers[ids,1]))>=3


def test_stencil_configuration_and_bad_widths(monkeypatch):
    assert environment({})['LAGFORCE_VELOCITY_STENCIL']=='area'
    assert environment({'velocity_stencil':'cell_width'})['LAGFORCE_VELOCITY_STENCIL']=='cell_width'
    with pytest.raises(ValueError):environment({'velocity_stencil':'bad'})
    monkeypatch.setenv('LAGFORCE_VELOCITY_STENCIL','bad')
    with pytest.raises(ValueError):stencil_mode()
    monkeypatch.setenv('LAGFORCE_VELOCITY_STENCIL','cell_width')
    with pytest.raises(ValueError):velocity_tensor_stencil(np.zeros((1,2)),np.zeros(2),1.,np.zeros((4,2)))

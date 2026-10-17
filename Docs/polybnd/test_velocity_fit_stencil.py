import os
import numpy as np
import pytest
from quadratic_velocity_reconstruction import curvature_fit,fit_mode
from pressure_fit_stencil import tensor_stencil
from physical_field_accuracy_audit import environment


def test_tensor_velocity_preserves_planar_fit_at_truncated_wall():
    x,y=np.meshgrid(np.linspace(-2,2,9),np.linspace(0,2,5));z=np.column_stack([x.ravel(),y.ravel()])
    ids,w=tensor_stencil(z,np.zeros(2),1.,3.)
    z=z[ids];A=np.column_stack([z[:,0]**2-.1,z[:,0]*z[:,1],z[:,1]**2-.1])
    target=np.column_stack([np.cos(z[:,0])-.4,np.zeros(len(z))])
    fit=curvature_fit(A,target,w,'tensor_free_mean')
    assert np.max(abs(fit[1:]))<1e-9
    assert np.max(abs(fit[:,1]))==0


def test_default_velocity_fit_matches_original_formula():
    rng=np.random.default_rng(17);A=rng.normal(size=(30,3));target=rng.normal(size=(30,2));w=rng.uniform(.1,1,30)
    expected=np.linalg.solve(A.T@(w[:,None]*A)+1e-10*np.eye(3),A.T@(w[:,None]*target))
    np.testing.assert_array_equal(curvature_fit(A,target,w,'radial_fixed_mean'),expected)


def test_velocity_fit_explicit_configuration(monkeypatch):
    monkeypatch.delenv('LAGFORCE_VELOCITY_FIT',raising=False)
    assert fit_mode()=='radial_fixed_mean'
    assert environment({})['LAGFORCE_VELOCITY_FIT']=='radial_fixed_mean'
    assert environment({'velocity_fit':'tensor_free_mean'})['LAGFORCE_VELOCITY_FIT']=='tensor_free_mean'
    monkeypatch.setenv('LAGFORCE_VELOCITY_FIT','bad')
    with pytest.raises(ValueError):fit_mode()
    with pytest.raises(ValueError):environment({'velocity_fit':'bad'})


def test_planar_odd_residual_needs_affine_nuisance_terms():
    x,y=np.meshgrid(np.linspace(-2,2,9),np.linspace(0,2,5));z=np.column_stack([x.ravel(),y.ravel()])
    ids,w=tensor_stencil(z,np.zeros(2),1.,3.);z=z[ids]
    A=np.column_stack([z[:,0]**2-.1,z[:,0]*z[:,1],z[:,1]**2-.1])
    target=np.column_stack([np.cos(z[:,0])+.3*z[:,0]-.4,np.zeros(len(z))])
    old=curvature_fit(A,target,w,'tensor_free_mean')
    assert abs(old[1,0])>.1
    fit=curvature_fit(A,target,w,'tensor_free_affine',zmean=z)
    assert np.max(abs(fit[1:]))<1e-9
    assert np.max(abs(fit[:,1]))==0
    assert environment({'velocity_fit':'tensor_free_affine'})['LAGFORCE_VELOCITY_FIT']=='tensor_free_affine'
    with pytest.raises(ValueError):curvature_fit(A,target,w,'tensor_free_affine')


@pytest.mark.parametrize('actual,expected',[
    ('radial_fixed_mean','velocity fit configuration mismatch'),
    ('tensor_free_affine','velocity fit requires quadratic reconstruction'),
])
def test_evolution_rejects_misreported_velocity_fit(monkeypatch,tmp_path,actual,expected):
    import json
    from simple_wave_evolution import run_with_executor
    monkeypatch.chdir(tmp_path)
    for key in ('LAGFORCE_PRESSURE_RIEMANN','LAGFORCE_VECTOR_KINETIC_ROOT','LAGFORCE_QUADRATIC_VELOCITY'):
        monkeypatch.delenv(key,raising=False)
    monkeypatch.setenv('LAGFORCE_VELOCITY_FIT',actual)
    (tmp_path/'simple_wave_config.json').write_text(json.dumps({'velocity_fit':'tensor_free_affine'}))
    with pytest.raises(ValueError,match=expected):run_with_executor()
    assert not (tmp_path/'initial_state.json').exists()
    assert not (tmp_path/'checkpoint.json').exists()

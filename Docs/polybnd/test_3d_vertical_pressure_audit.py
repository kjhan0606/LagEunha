import numpy as np
from audit_3d_vertical_pressure import gradient_z
from pressure_reconstruction_3d import value


def test_gradient_matches_polynomial_derivative():
    rng = np.random.default_rng(317)
    metadata = dict(center=rng.normal(size=3), scale=.7, mean=2.,
                    coefficient=rng.normal(size=9), feature_mean=rng.normal(size=9))
    points = rng.normal(size=(30, 3))
    shift = np.array([0., 0., 1e-5])
    numerical = (value(metadata, points+shift)-value(metadata, points-shift))/2e-5
    np.testing.assert_allclose(gradient_z(metadata, points), numerical, atol=1e-8, rtol=1e-8)


def test_planar_polynomial_has_zero_vertical_gradient():
    metadata = dict(center=np.zeros(3), scale=.2,
                    coefficient=np.array([1., 2., 0., 3., 4., 0., 5., 0., 0.]))
    np.testing.assert_array_equal(gradient_z(metadata, np.ones((4, 3))), np.zeros(4))

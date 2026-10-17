import unittest
import numpy as np
from test_subcell_moments import setup_case
from subcell_moments import transfer, moments
from subcell_quadrature import quadrature, features, compress, subdivide


class QuadratureTests(unittest.TestCase):
    def test_positive_bounded_moment_representation(self):
        for kind in ('rotation', 'shear', 'kepler'):
            for ratio in (1., .01):
                ff, _, receivers = setup_case(kind, ratio)
                ff = subdivide(ff)
                groups, _ = transfer(ff, receivers)
                # Include combined patch so the elimination is exercised too.
                for group in groups+[ff]:
                    points, weights = quadrature(group)
                    target = features(points)@weights
                    p, w = compress(points, weights)
                    self.assertLessEqual(len(w), 16)
                    self.assertTrue((w > 0).all())
                    np.testing.assert_allclose(features(p)@w, target, rtol=3e-12, atol=2e-15)
                    m = moments(group)
                    self.assertAlmostEqual(target[9]-target[10], m['angular'], places=14)
                    self.assertAlmostEqual((target[12]+target[14])/2, m['kinetic'], places=14)
                    self.assertAlmostEqual(target[15], m['internal'], places=15)

    def test_velocity_jump_keeps_internal_energy_separate(self):
        ff, _, receivers = setup_case('kepler', .01)
        for f in ff[3:]: f.b = f.b+np.array([.7, -.2])
        groups, _ = transfer(ff, receivers)
        for group in groups:
            p, w = quadrature(group); target = features(p)@w
            p, w = compress(p, w)
            np.testing.assert_allclose(features(p)@w, target, rtol=3e-12, atol=2e-15)
            self.assertTrue((p[:, -1] > 0).all())

    def test_fixed_kepler_potential_additional_moment(self):
        ff, _, receivers = setup_case('kepler', .01)
        ff = subdivide(ff); groups, _ = transfer(ff, receivers)
        potential = lambda x: -1/np.linalg.norm(x, axis=1)
        for group in groups+[ff]:
            p, w = quadrature(group); target = features(p, potential)@w
            p, w = compress(p, w, potential)
            self.assertLessEqual(len(w), 17)
            self.assertTrue((w > 0).all())
            np.testing.assert_allclose(features(p, potential)@w, target, rtol=3e-12, atol=2e-15)


if __name__ == '__main__':
    unittest.main(verbosity=2)

import unittest
import numpy as np
from mass_feasibility import minimum_kinetic, audit


class Feasibility(unittest.TestCase):
    def test_rigid_translation_rotation_bound(self):
        x = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
        # M=4, centroid=(.5,.5), I=2, P=(4,0), Lcentroid=2.
        # Kmin=P^2/(2M)+Lcentroid^2/(2I)=2+1=3.
        self.assertAlmostEqual(minimum_kinetic(x, np.ones(4), np.array([4., 0.]), 0.), 3.)

    def test_cocircular_cold_contact_counterexample(self):
        for reverse in (False, True):
            result = audit(.1, .01, 2.78e-5, 0., reverse)
            self.assertFalse(result['feasible_positive_internal'])
            self.assertLess(result['maximum_possible_internal'], -5e-5)
            self.assertLess(result['maximum_internal_excluding_gravity'], 0.)

    def test_uniform_density_is_not_blanket_rejected(self):
        result = audit(.1, 1., 2.78e-5, 0.)
        self.assertTrue(result['feasible_positive_internal'])


if __name__ == '__main__':
    unittest.main(verbosity=2)

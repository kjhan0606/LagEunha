import unittest
import numpy as np
from fixed_mass_allocation import audit


class FixedMassTests(unittest.TestCase):
    def test_square_flip_has_fixed_mass_allocation(self):
        square = [[-.1, -.1], [.1, -.1], [.1, .1], [-.1, .1]]
        for rho in ((1., 1.), (1., .01), (.01, 1.)):
            result = audit(square, rho)
            self.assertTrue(result['feasible'])
            self.assertLess(result['maximum_constraint_residual'], 1e-12)

    def test_cocircular_irregular_flip_can_be_infeasible(self):
        # All four vertices lie on the unit circle; 24^2+7^2=25^2.
        x = np.array([[1., 0.], [24/25, 7/25], [-1., 0.], [0., -1.]])
        np.testing.assert_allclose(np.sum(x*x, axis=1), 1., atol=1e-15)
        result = audit(x, (1., 1.))
        self.assertFalse(result['feasible'])
        self.assertGreater(result['degree_one_deficit'], .26)


if __name__ == '__main__':
    unittest.main(verbosity=2)

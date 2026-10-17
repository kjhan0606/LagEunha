import unittest
import numpy as np
from feedback_patch import initial, mesh, geom, pressure_force, cross, GAMMA, run


class Conservation(unittest.TestCase):
    def test_force_is_energy_gradient(self):
        x = initial(.1); keys = mesh(x); mass, _ = geom(x, keys)
        p0 = .278; entropy = np.full(len(keys), p0)
        x[0] += [.003, -.002]
        force, _, _ = pressure_force(x, keys, mass, entropy, p0)
        def potential(x):
            area, _ = geom(x, keys)
            return (entropy*mass**GAMMA*area**(1-GAMMA)/(GAMMA-1)+p0*area).sum()
        numerical = np.zeros_like(x); eps = 1e-7
        for i in range(4):
            for j in range(2):
                dx = np.zeros_like(x); dx[i, j] = eps
                numerical[i, j] = -(potential(x+dx)-potential(x-dx))/(2*eps)
        np.testing.assert_allclose(force, numerical, rtol=2e-7, atol=1e-10)
        np.testing.assert_allclose(force.sum(0), 0., atol=1e-15)
        self.assertLess(abs(cross(x, force).sum()), 1e-15)

    def test_uniform_pressure_reservoir_balance(self):
        x = initial(.1); keys = mesh(x); mass, _ = geom(x, keys)
        force, _, _ = pressure_force(x, keys, mass, np.full(len(keys), .278), .278)
        np.testing.assert_allclose(force, 0., atol=1e-15)

    def test_closed_energy_and_disabled_force_work(self):
        for feedback in (0, 1):
            result, _, _ = run(.1, .278, .001, feedback, False, duration=.04)
            self.assertLess(abs(result['energy_residual']), 1e-10)
            self.assertLess(abs(result['lz_relative']), 1e-10)
            if feedback == 0:
                self.assertGreater(abs(result['control_work']), 1e-10)


if __name__ == '__main__':
    unittest.main(verbosity=2)

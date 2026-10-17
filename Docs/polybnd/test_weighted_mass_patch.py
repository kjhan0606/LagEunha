import unittest
import numpy as np
from weighted_mass_patch import CapacityPatch, extended_support


SQUARE = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])
IRREGULAR = np.array([[1., 0.], [24/25, 7/25], [-1., 0.], [0., -1.]])


class WeightTests(unittest.TestCase):
    def test_extended_support_preserves_fixed_labels(self):
        for x in (SQUARE, IRREGULAR):
            for rho in ((1., 1.), (1., .01), (.01, 1.)):
                result = extended_support(x, rho)
                patch = CapacityPatch(x, rho)
                np.testing.assert_allclose(result['allocation'].sum(1), patch.fixed, rtol=1e-13, atol=1e-15)
                self.assertAlmostEqual(result['allocation'].sum(), patch.total, places=13)
        result = extended_support(IRREGULAR, (1., 1.))
        self.assertGreater(result['nonincident'], .266)

    def test_weights_solve_original_infeasible_patch(self):
        for x in (SQUARE, IRREGULAR):
            for rho in ((1., 1.), (1., .01), (.01, 1.)):
                patch = CapacityPatch(x, rho); sites = x.copy()
                weights, _, error = patch.solve(sites)
                np.testing.assert_array_equal(sites, x)
                mass, volume, cells = patch.evaluate(sites, weights)
                np.testing.assert_allclose(mass, patch.fixed, rtol=1e-10, atol=1e-13)
                self.assertLess(error, 1e-10)
                self.assertTrue((volume > 0).all())
                self.assertAlmostEqual(volume.sum(), sum(p.area for p in patch.donors), places=12)

    def test_weight_gauge(self):
        patch = CapacityPatch(IRREGULAR)
        w, _, _ = patch.solve(IRREGULAR)
        np.testing.assert_allclose(patch.evaluate(IRREGULAR, w)[0],
                                   patch.evaluate(IRREGULAR, w+3.2)[0], atol=1e-13)

    def test_weight_changes_for_changed_sites_with_fixed_mass(self):
        patch = CapacityPatch(SQUARE, (1., .01)); results = []
        direction = np.array([[.3, .1], [-.1, .2], [-.2, -.1], [.2, -.2]])
        for parameter in (-.1, 0., .1):
            sites = .8*SQUARE+parameter*direction
            w, _, _ = patch.solve(sites); results.append(w)
            np.testing.assert_allclose(patch.evaluate(sites, w)[0], patch.fixed, rtol=1e-10)
        self.assertGreater(np.linalg.norm(results[0]-results[-1]), 1e-3)

    def test_pressure_work_includes_weights(self):
        patch = CapacityPatch(SQUARE)
        sites = .8*SQUARE; w, _, _ = patch.solve(sites)
        mass, area, _ = patch.evaluate(sites, w)
        gamma = 5/3; p = np.array([1., 2., .5, 1.3]); entropy = p/(mass/area)**gamma
        def U(weights):
            a = patch.evaluate(sites, weights)[1]
            return np.sum(entropy*mass**gamma*a**(1-gamma)/(gamma-1))
        jac = patch.volume_jacobian(sites, w)
        direction = np.array([.3, -.2, .1, -.1]); eps = 1e-6
        fd = (U(w+eps*direction)-U(w-eps*direction))/(2*eps)
        analytic = -p@jac@direction
        self.assertGreater(abs(fd), .01)
        self.assertAlmostEqual(fd, analytic, places=7)
        np.testing.assert_allclose(jac.sum(0), 0., atol=1e-9)
        np.testing.assert_allclose(jac, jac.T, atol=2e-6)

    def test_capacity_constraint_chain_rule(self):
        patch = CapacityPatch(SQUARE); sites = .8*SQUARE
        w, _, _ = patch.solve(sites)
        direction = np.array([[.3, .1], [-.1, .2], [-.2, -.1], [.2, -.2]])
        eps = 2e-5; plus = sites+eps*direction; minus = sites-eps*direction
        wp, _, _ = patch.solve(plus); wm, _, _ = patch.solve(minus)
        fixed_w_derivative = (patch.evaluate(plus, w)[1]-patch.evaluate(minus, w)[1])/(2*eps)
        weight_derivative = patch.volume_jacobian(sites, w)@((wp-wm)/(2*eps))
        np.testing.assert_allclose(fixed_w_derivative+weight_derivative, 0., atol=2e-8)
        pressure = np.array([1., 2., .5, 1.3])
        self.assertGreater(abs(pressure@fixed_w_derivative), .01)

    def test_zero_mass_flux_does_not_remove_momentum_energy_flux(self):
        # Rigid rotation v=(-y,x). A square with centroid (1,0) translates
        # at the material velocity of its generator (.75,-.1): no steering.
        generator = np.array([.75, -.1]); face_velocity = np.array([.1, .75])
        nodes, weights = np.polynomial.legendre.leggauss(2)
        mass_flux = 0.; momentum_flux = np.zeros(2); kinetic_flux = 0.
        faces = [(np.array([.5, 0.]), np.array([0., .5]), np.array([-1., 0.])),
                 (np.array([1.5, 0.]), np.array([0., .5]), np.array([1., 0.])),
                 (np.array([1., -.5]), np.array([.5, 0.]), np.array([0., -1.])),
                 (np.array([1., .5]), np.array([.5, 0.]), np.array([0., 1.]))]
        for center, tangent, normal in faces:
            for q, weight in zip(nodes, weights):
                x = center+q*tangent; v = np.array([-x[1], x[0]])
                flux = np.dot(v-face_velocity, normal)*weight*np.linalg.norm(tangent)
                mass_flux += flux; momentum_flux += flux*v
                kinetic_flux += flux*np.dot(v, v)/2
        self.assertAlmostEqual(mass_flux, 0., places=14)
        np.testing.assert_allclose(momentum_flux, [-.25, -.1], atol=1e-14)
        self.assertAlmostEqual(kinetic_flux, -.1, places=14)

    def test_actual_power_cell_ale_balance(self):
        from weighted_ale_balance import audit
        result = audit()
        self.assertLess(result['max_balance_error'], 3e-7)
        self.assertLess(result['mass_flux_max'], 3e-7)
        self.assertGreater(result['advective_momentum_max'], .01)
        self.assertGreater(result['omit_advective_error'], .01)
        self.assertGreater(result['weight_rate_norm'], .001)


if __name__ == '__main__':
    unittest.main(verbosity=2)

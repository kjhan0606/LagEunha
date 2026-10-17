import unittest
import numpy as np
from enriched_flip import EnrichedFlip, cross


def kepler_velocity(x):
    return np.column_stack((-x[:, 1], x[:, 0]))/np.linalg.norm(x, axis=1)[:, None]**1.5


class EnrichedTests(unittest.TestCase):
    def test_single_flip_preserves_nodal_velocity_and_invariants(self):
        for rho in ((1., 1.), (1., .01), (.01, 1.)):
            p = EnrichedFlip(rho=rho); v = kepler_velocity(p.x); V = p.encode_velocity(v)
            s = p.state(p.Q0, V); oldmom = p.Mold@v
            np.testing.assert_allclose(s['physical_velocity'][:4], v, atol=0.)
            np.testing.assert_allclose(s['momentum'], oldmom.sum(0), rtol=1e-13, atol=1e-15)
            self.assertAlmostEqual(s['angular'], np.sum(cross(p.x, oldmom)), places=14)
            self.assertAlmostEqual(s['kinetic'], .5*np.sum(v*oldmom), places=14)
            self.assertGreater(s['internal'], 0.)
            self.assertGreater(np.linalg.eigvalsh(p.M).min(), 0.)

    def test_affine_flow_needs_no_internal_velocity(self):
        p = EnrichedFlip()
        v = p.x@np.array([[.2, -1.7], [1.7, -.3]]).T+np.array([.5, -.2])
        np.testing.assert_allclose(p.encode_velocity(v)[4], 0., atol=1e-15)

    def test_cell_interior_bubble_is_insufficient(self):
        p = EnrichedFlip(); v = kepler_velocity(p.x); V = p.encode_velocity(v)
        center = (p.T@p.Q0)[4]
        mismatch = p.residual_value(center, p.Q0, v)
        self.assertGreater(np.linalg.norm(mismatch), 1e-4)
        np.testing.assert_allclose(mismatch, V[4], rtol=1e-12, atol=1e-15)
        # A bubble zero on every new triangle edge is zero on the new diagonal.
        # It cannot represent this nonzero mismatch at the crossing.

    def test_old_field_reproduced_everywhere_by_crossing_mode(self):
        p = EnrichedFlip(); v = kepler_velocity(p.x); V = p.encode_velocity(v)
        y = p.T@p.Q0; vel = p.T@V
        for ids in p.tri:
            for bary in ([.2, .3, .5], [.6, .3, .1]):
                point = np.array(bary)@y[list(ids)]
                represented = np.array(bary)@vel[list(ids)]
                for oldids in p.old_tri:
                    lam = np.linalg.solve(np.vstack((np.ones(3), p.x[list(oldids)].T)), np.r_[1., point])
                    if lam.min() >= -1e-12:
                        np.testing.assert_allclose(represented, lam@v[list(oldids)], atol=1e-14)
                        break
                else: self.fail('no old triangle found')

    def test_variational_force_work_and_symmetries(self):
        p = EnrichedFlip(pressure=.278); Q = p.Q0.copy(); Q[4] += [.004, -.003]
        Q[0] += [-.002, .001]
        force = p.force(Q); eps = 1e-7; fd = np.zeros_like(Q)
        for i in range(5):
            for j in range(2):
                dq = np.zeros_like(Q); dq[i, j] = eps
                fd[i, j] = -(p.potential(Q+dq)-p.potential(Q-dq))/(2*eps)
        np.testing.assert_allclose(force, fd, rtol=3e-7, atol=1e-10)
        np.testing.assert_allclose(force[:4].sum(0), 0., atol=1e-15)
        self.assertLess(abs(np.sum(cross(Q, force))), 1e-15)
        V = np.random.default_rng(123).normal(size=(5, 2))
        udot = (p.potential(Q+eps*V)-p.potential(Q-eps*V))/(2*eps)
        kdot = np.sum(V*(p.M@p.acceleration(Q)))
        self.assertLess(abs(udot+kdot), 1e-9)

    def test_inertia_cross_terms_are_required(self):
        p = EnrichedFlip(); self.assertGreater(np.linalg.norm(p.M[:4, 4]), 1e-5)
        v = kepler_velocity(p.x); V = p.encode_velocity(v)
        exact = .5*np.sum(V*(p.M@V))
        split = .5*np.sum(V[:4]*(p.M[:4, :4]@V[:4]))+.5*p.M[4, 4]*np.dot(V[4], V[4])
        self.assertGreater(abs(exact-split), 1e-7)


if __name__ == '__main__':
    unittest.main(verbosity=2)

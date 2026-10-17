import unittest
import numpy as np
from subcell_exact_geometry import Polygon
from subcell_moments import (Fragment, GAMMA, polygon_moments, corner_polygons,
                             transfer, moments, invariants)


def setup_case(kind='rotation', contrast=1.):
    x = np.array([[.7, -.1], [.9, -.1], [.9, .1], [.7, .1]])
    old = [(0, 1, 2), (0, 2, 3)]; new = [(0, 1, 3), (1, 2, 3)]
    donors, labels = corner_polygons(x, old)
    receivers, _ = corner_polygons(x, new)
    fields = []
    for ids in old:
        if kind == 'rotation':
            B = np.array([[0., -1.7], [1.7, 0.]]); b = np.array([.3, -.8])
        elif kind == 'shear':
            B = np.array([[0., 2.], [0., 0.]]); b = np.array([.4, -.1])
        else:
            p = x[list(ids)]; r = np.linalg.norm(p, axis=1)
            v = np.column_stack((-p[:, 1], p[:, 0]))/r[:, None]**1.5
            coeff = np.linalg.solve(np.column_stack((np.ones(3), p)), v)
            b = coeff[0]; B = coeff[1:].T
        fields.append((b, B))
    fragments = [Fragment(p, 1. if t == 0 else contrast, 2.78e-5/(GAMMA-1),
                          fields[t][0], fields[t][1]) for p, (t, _) in zip(donors, labels)]
    return fragments, donors, receivers


class SubcellTests(unittest.TestCase):
    def test_polygon_integrals(self):
        area, q, Q = polygon_moments(Polygon([(0, 0), (1, 0), (0, 1)]))
        self.assertAlmostEqual(area, .5)
        np.testing.assert_allclose(q, [1/6, 1/6], atol=1e-15)
        np.testing.assert_allclose(Q, [[1/12, 1/24], [1/24, 1/12]], atol=1e-15)

    def test_affine_rotation_and_shear(self):
        for kind in ('rotation', 'shear'):
            for contrast in (1., .01):
                fragments, _, receivers = setup_case(kind, contrast)
                groups, closure = transfer(fragments, receivers)
                after = [f for g in groups for f in g]
                np.testing.assert_allclose(invariants(after), invariants(fragments), rtol=2e-13, atol=2e-15)
                self.assertLess(closure, 1e-12)
                for g in groups:
                    m = moments(g)
                    np.testing.assert_allclose(m['B'], fragments[0].B, rtol=1e-11, atol=1e-11)
                    self.assertLess(abs(m['residual_energy']), 1e-15)
                    self.assertAlmostEqual(m['pressure'], 2.78e-5, places=16)

    def test_corner_spin_is_required(self):
        fragments, _, receivers = setup_case()
        groups, _ = transfer(fragments, receivers)
        total = moments(fragments); orbital = 0.; spin = 0.
        for g in groups:
            m = moments(g); c = m['center']; p = m['momentum']
            orbital += c[0]*p[1]-c[1]*p[0]; spin += m['spin']
        self.assertGreater(abs(spin), 1e-6)
        self.assertAlmostEqual(orbital+spin, total['angular'], places=14)

    def test_repeated_flip_kepler_contact(self):
        fragments, old, new = setup_case('kepler', .01)
        baseline = invariants(fragments)
        for _ in range(20):
            groups, _ = transfer(fragments, new)
            for g in groups:
                m = moments(g)
                self.assertGreater(m['internal'], 0.)
                self.assertGreater(m['mass'], 0.)
                self.assertGreater(m['residual_min_eigenvalue'], -1e-14)
                self.assertAlmostEqual(m['pressure'], 2.78e-5, places=15)
            fragments = [f for g in groups for f in g]; old, new = new, old
        np.testing.assert_allclose(invariants(fragments), baseline, rtol=1e-12, atol=1e-15)

    def test_velocity_jump_retains_nonthermal_energy(self):
        fragments, _, receivers = setup_case('shear', .01)
        for i, f in enumerate(fragments):
            if i >= 3: f.b = f.b+np.array([.7, -.2])
        groups, _ = transfer(fragments, receivers)
        residual = sum(moments(g)['residual_energy'] for g in groups)
        self.assertGreater(residual, 1e-8)
        np.testing.assert_allclose(invariants([f for g in groups for f in g]), invariants(fragments), rtol=1e-12, atol=1e-15)
        for g in groups: self.assertAlmostEqual(moments(g)['pressure'], 2.78e-5, places=15)

    def test_quadratic_potential_integral(self):
        fragments, _, receivers = setup_case('kepler', .01)
        groups, _ = transfer(fragments, receivers)
        H = np.array([[2., .3], [.3, -1.]])
        def potential(ff):
            m = moments(ff); Q = m['S']+m['mass']*np.outer(m['center'], m['center'])
            return 3*m['mass']+np.dot([1., -2.], m['first'])+np.sum(H*Q)/2
        self.assertAlmostEqual(sum(potential(g) for g in groups), potential(fragments), places=14)

    def test_missing_and_duplicate_receivers_rejected(self):
        fragments, _, receivers = setup_case()
        with self.assertRaisesRegex(RuntimeError, 'coverage'):
            transfer(fragments, receivers[:-1])
        with self.assertRaisesRegex(RuntimeError, 'overlapping'):
            transfer(fragments, receivers+[receivers[0]])


if __name__ == '__main__':
    unittest.main(verbosity=2)

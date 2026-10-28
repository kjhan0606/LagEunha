"""Independent geometric rejection examples for the W=0 certificate."""
import unittest

import numpy as np

from audit_kepler_vertex_halo_certificate_20261029 import influence_bounds, rectangle_margin


class VertexCertificateTests(unittest.TestCase):
    def test_elongated_exact_cell_passes_tighter_rectangle(self):
        owner = np.zeros(2)
        vertices = np.array([[-.1, -1.], [.1, -1.], [.1, 1.], [-.1, 1.]])
        lower, upper = np.array([-1.5, -2.1]), np.array([1.5, 2.1])
        lo, hi, radius = influence_bounds(owner, vertices)
        self.assertGreater(rectangle_margin((lo, hi), lower, upper), 0.)
        self.assertLess(rectangle_margin((owner-2*radius, owner+2*radius), lower, upper), 0.)
        # This rectangle is the exact Voronoi cell for these four sites.
        sites = np.array([[-.2, 0.], [.2, 0.], [0., -2.], [0., 2.]])
        gaps = np.sum((vertices[:, None] - sites)**2, axis=2) - np.sum(vertices**2, axis=1)[:, None]
        self.assertTrue(np.all(gaps >= -1e-15))
        # Directly check bisector inequalities for a grid of omitted sites;
        # no use of the certificate formula in this independent oracle.
        grid = np.linspace(-4., 4., 161)
        x, y = np.meshgrid(grid, grid)
        candidates = np.column_stack((x.ravel(), y.ravel()))
        omitted = candidates[np.any((candidates < lower) | (candidates > upper), axis=1)]
        gaps = np.sum((vertices[:, None] - omitted)**2, axis=2) - np.sum(vertices**2, axis=1)[:, None]
        self.assertGreaterEqual(float(np.min(gaps)), 0.)

    def test_insufficient_rectangle_rejected_and_real_cutter_exists(self):
        vertices = np.array([[-.1, -1.], [.1, -1.], [.1, 1.], [-.1, 1.]])
        lower, upper = np.array([-.5, -2.1]), np.array([.5, 2.1])
        bounds = influence_bounds(np.zeros(2), vertices)
        self.assertLess(rectangle_margin(bounds, lower, upper), 0.)
        omitted_site = np.array([.6, 1.])
        gaps = np.sum((vertices - omitted_site)**2, axis=1) - np.sum(vertices**2, axis=1)
        self.assertLess(float(np.min(gaps)), 0.)

    def test_coordinate_translation_preserves_certificate(self):
        owner = np.array([.3, .2])
        vertices = np.array([[-.1, -1.], [.1, -1.], [.1, 1.], [-.1, 1.]]) + owner
        lower, upper = np.array([-1.5, -2.1]) + owner, np.array([1.5, 2.1]) + owner
        original = rectangle_margin(influence_bounds(owner, vertices), lower, upper)
        shift = np.array([6., -6.])
        translated = rectangle_margin(influence_bounds(owner+shift, vertices+shift), lower+shift, upper+shift)
        self.assertAlmostEqual(original, translated, places=14)


if __name__ == '__main__':
    unittest.main()

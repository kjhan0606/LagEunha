import unittest
import numpy as np
from moment_aligned_density_audit import polygon_moments
from density_moment_bounds import density_moment_bounds


class DensityMomentBoundTests(unittest.TestCase):
    def geometry(self):return polygon_moments(np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]]))

    def test_constant_density_and_translation(self):
        for shift in (np.zeros(2),np.array([3.,-2.])):
            vertices=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])+shift;g=polygon_moments(vertices)
            bounds=density_moment_bounds(g,4.2*g)
            self.assertAlmostEqual(bounds['density_eigen_min'],4.2,places=11)
            self.assertAlmostEqual(bounds['density_eigen_max'],4.2,places=11)

    def test_half_contact_requires_more_than_cell_mean_upper_bound(self):
        g=self.geometry();high=polygon_moments(np.array([[.5,0.],[1.,0.],[1.,1.],[.5,1.]]))
        target=g+99*high;b=density_moment_bounds(g,target)
        self.assertAlmostEqual(b['density_eigen_max'],50.5+99*np.sqrt(3)/4,places=11)
        self.assertAlmostEqual(b['density_eigen_min'],50.5-99*np.sqrt(3)/4,places=11)
        self.assertGreater(b['density_eigen_max'],1.01*target[0]/g[0])

    def test_negative_input_is_diagnosed_without_projection(self):
        g=self.geometry();b=density_moment_bounds(g,-g)
        self.assertAlmostEqual(b['density_eigen_max'],-1.)
        with self.assertRaises(ValueError):density_moment_bounds(np.zeros(7),g)


if __name__=='__main__':unittest.main()

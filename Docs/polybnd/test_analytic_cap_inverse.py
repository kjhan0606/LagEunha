import unittest
import numpy as np
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon
from analytic_cap_inverse import cap_at_fraction
from moment_aligned_density_audit import cut_at_fraction


class AnalyticCapTests(unittest.TestCase):
    def test_experimental_dispatch_is_explicit(self):
        import os
        from unittest.mock import patch
        vertices=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        normal=np.array([1.,.3]);expected_offset,expected=cap_at_fraction(vertices,normal,.23)
        with patch.dict(os.environ,{'LAGFORCE_ANALYTIC_CAP':'1'}):
            offset,actual=cut_at_fraction(vertices,normal,.23)
        self.assertEqual(offset,expected_offset)
        np.testing.assert_array_equal(actual.exterior.coords,expected.exterior.coords)

    def test_square_vertex_ties_and_complements(self):
        vertices=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        for normal in (np.array([1.,0.]),np.array([0.,1.]),np.array([1.,1.])):
            for fraction in (0.,1e-12,1e-7,.1,.5,.9,1-1e-7,1-1e-12,1.):
                _,cut=cap_at_fraction(vertices,normal,fraction)
                self.assertLess(abs(cut.area-fraction),2e-13)

    def test_random_convex_caps_match_target_and_reference(self):
        rng=np.random.default_rng(20261016)
        for shift in (np.zeros(2),np.array([3.,-2.])):
            for _ in range(8):
                points=rng.uniform(-1,1,(12,2));vertices=points[ConvexHull(points).vertices]+shift
                area=Polygon(vertices).area
                for theta in rng.uniform(-np.pi,np.pi,4):
                    normal=np.array([np.cos(theta),np.sin(theta)])
                    for fraction in (1e-8,.03,.3,.7,.97,1-1e-8):
                        _,actual=cap_at_fraction(vertices,normal,fraction)
                        _,reference=cut_at_fraction(vertices,normal,fraction)
                        self.assertLess(abs(actual.area/area-fraction),2e-13)
                        self.assertLess(actual.symmetric_difference(reference).area/area,3e-12)

    def test_invalid_fraction_and_normal(self):
        square=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        for fraction in (-.1,1.1,float('nan')):
            with self.assertRaises(ValueError):cap_at_fraction(square,np.ones(2),fraction)
        with self.assertRaises(ValueError):cap_at_fraction(square,np.zeros(2),.5)


if __name__=='__main__':unittest.main()

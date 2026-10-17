import unittest
import numpy as np
from shapely.geometry import Polygon
from shear_interface_audit import moments
from subcell_density_audit import triangles,triangle_moments,fit


class SubcellTests(unittest.TestCase):
    def test_partition_moments(self):
        p=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        for level in (1,2):
            total=sum((triangle_moments(t) for t in triangles(p,level)),np.zeros(7))
            np.testing.assert_allclose(total,moments(Polygon(p)),atol=1e-14,rtol=0)

    def test_constant_feasible(self):
        p=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        result,_,rho=fit(p,2*moments(Polygon(p)),1,lambda x:np.full(len(x),2.))
        self.assertTrue(result['moment_gate']);self.assertLess(result['mass_relative_error'],1e-10)
        self.assertTrue(np.all(rho>=0))

    def test_outside_centroid_rejected(self):
        p=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        target=moments(Polygon(p));target[1]=2.
        result,_,_=fit(p,target,1,lambda x:np.ones(len(x)))
        self.assertFalse(result['feasible'])

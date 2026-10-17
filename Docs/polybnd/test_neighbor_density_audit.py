import unittest
import numpy as np
from neighbor_density_audit import interface_prior
from discontinuous_density_audit import clip


class NeighborTests(unittest.TestCase):
    def test_normal_and_mass_fraction(self):
        vertices=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        centers=np.array([[.5,.5],[-.5,.5],[1.5,.5],[.5,-.5],[.5,1.5]])
        means=np.array([5.5,1,10,5.5,5.5]);neighbors=[{1,2,3,4},set(),set(),set(),set()]
        seed,d=interface_prior(0,vertices,centers,means,neighbors)
        np.testing.assert_allclose(d['normal'],[1,0],atol=1e-14)
        self.assertAlmostEqual(clip(vertices,np.array(d['normal']),d['offset']).area,.5,places=12)
        np.testing.assert_array_equal(seed(np.array([[.1,.5],[.9,.5]])),[1,10])

    def test_uniform_field(self):
        p=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        seed,d=interface_prior(0,p,np.array([[.5,.5]]),np.array([2.]),[set()])
        self.assertTrue(d['flat']);np.testing.assert_array_equal(seed(p),2.)

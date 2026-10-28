"""Small geometric tests; no simulation or large tessellation."""
import unittest
import numpy as np
from audit_kepler_t3_support_20261029 import geometry, orbit


class GeometryTests(unittest.TestCase):
    def setUp(self):
        self.xy=np.array([(x,y) for x in range(1,6) for y in range(1,6)],dtype=float)
        self.selected=np.array([12])
        self.foreign=np.zeros(25,dtype=bool)
        self.foreign[17]=True

    def test_periodic_certificate_matches_explicit_images(self):
        area, support, pairs, safety, _=geometry(self.xy,self.selected,self.foreign,True)
        self.assertAlmostEqual(area[0],1.)
        self.assertTrue(support[0])
        self.assertIn((12,17),pairs)
        self.assertGreater(safety,1.)
        copies=np.vstack([self.xy]+[self.xy+6*np.array([i,j])
            for i in [-1,0,1] for j in [-1,0,1] if (i,j)!=(0,0)])
        explicit,_,_,_,_=geometry(copies,self.selected,np.tile(self.foreign,9),False)
        np.testing.assert_allclose(area,explicit,atol=1.e-14,rtol=0.)

    def test_outside_sites_rejected_for_periodic_certificate(self):
        self.xy[0]=[-.1,1.]
        with self.assertRaisesRegex(ValueError,'every generator'):
            geometry(self.xy,self.selected,self.foreign,True)

    def test_uncertified_edge_polygon_rejected(self):
        with self.assertRaises(ValueError):
            geometry(self.xy,np.array([0]),self.foreign,True)

    def test_exact_orbit_preserves_radius(self):
        evolved=orbit(self.xy,3.)
        np.testing.assert_allclose(np.linalg.norm(evolved-3.,axis=1),
                                   np.linalg.norm(self.xy-3.,axis=1),atol=1.e-14,rtol=0.)


if __name__=='__main__':unittest.main()

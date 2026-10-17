import unittest
import numpy as np
from shapely.geometry import Polygon
from discontinuous_density_audit import clip,fit
from shear_interface_audit import moments


class DiscontinuousTests(unittest.TestCase):
    def test_clip_and_recover(self):
        p=Polygon([(0,0),(1,0),(1,1),(0,1)]);v=np.asarray(p.exterior.coords[:-1])
        cut=clip(v,np.array([1.,0.]),.4)
        self.assertAlmostEqual(cut.area,.6)
        target=moments(p)+9*moments(cut)
        result,recovered=fit(p,target,1.,10.)
        self.assertTrue(result['feasible']);self.assertTrue(result['moment_gate'])
        self.assertLess(result['mass_relative_error'],1e-12)
        self.assertLess(cut.symmetric_difference(recovered).area,1e-7)

    def test_impossible_mass(self):
        p=Polygon([(0,0),(1,0),(1,1),(0,1)])
        result,_=fit(p,20*moments(p),1.,10.)
        self.assertFalse(result['feasible'])

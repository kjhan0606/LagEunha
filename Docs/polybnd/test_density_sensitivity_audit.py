import unittest
import numpy as np
from shapely.geometry import Polygon
from density_sensitivity_audit import rotated_prior
from discontinuous_density_audit import clip


class RotationTests(unittest.TestCase):
    def test_rotation_keeps_mass_fraction(self):
        v=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        metadata=dict(flat=False,normal=[1.,0.],offset=.7,fraction=.3)
        for angle in (-15.,15.):
            _,m=rotated_prior(v,None,metadata,1.,100.,angle)
            self.assertAlmostEqual(clip(v,np.array(m['normal']),m['offset']).area/Polygon(v).area,.3,places=12)
            theta=np.degrees(np.arctan2(m['normal'][1],m['normal'][0]))
            self.assertAlmostEqual(theta,angle,places=12)

    def test_flat_prior_is_preserved(self):
        prior=lambda x:np.ones(len(x))
        result,meta=rotated_prior(None,prior,dict(flat=True),1.,100.,15.)
        self.assertIs(result,prior);self.assertTrue(meta['flat'])


if __name__=='__main__':unittest.main()

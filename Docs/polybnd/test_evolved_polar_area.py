import unittest
import numpy as np
from scipy.spatial import ConvexHull
from kepler_disk_reference import polar_rule


class EvolvedPolarAreaTests(unittest.TestCase):
    def test_evolved_near_origin_face_layer(self):
        p=np.array(POLYGON)
        expected=ConvexHull(p-p[0]).volume
        for order in (16,24,32):
            x,w=polar_rule(p,order,radii=(.25,.5,2.,2.5))
            self.assertLess(abs(w.sum()-expected)/expected,2e-10)
            self.assertTrue(np.all(w>0))

POLYGON=[[3.184783199232708, 2.8072608377566164], [3.2827628103870334, 2.90319117363017], [3.2625317777478378, 3.0423880943206743], [3.2530512677431807, 3.050134619474138], [3.0000352859196915, 2.999960460016231], [3.169268475214074, 2.8067632270781124]]

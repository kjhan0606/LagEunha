import unittest
from unittest.mock import patch
import numpy as np
from potential_enriched_density import PotentialDensity
from kepler_disk_reference import polar_rule
from audit_kepler_energy_step import phi

FIXTURE={'cell': 13, 'polygon': [[3.4598873276045543, 3.145466161918769], [3.3695223676663946, 3.2730788520937453], [3.3684351799316152, 3.273269186227579], [3.23946143665229, 3.1772306943775144], [3.2530512677431807, 3.0501346194741386], [3.2625317777478378, 3.0423880943206743], [3.436061632986542, 3.0697526724370294]], 'raw': [0.0003444321228227567, 0.0011521370861150434, 0.0010831743945761737, 0.0038549348987928843, 0.003623372982687523, 0.003623372982687523, 0.0034073717840691093, -0.00020726213616110715, 0.0004967357296776962, -0.0006919761510605084, 0.0016611562336642943, -0.0006549133350118597, 0.0015603157976924267, 0.0001379887870603421, -0.0002935305630480075, -0.0002935305630480075, 0.0007199939275151391, 0.00043041140317858], 'coefficient': [-2.6686064130224723, 1.1563267049891488, 0.4783825969688036, -0.35433330722106915, -0.5271941017674795, 0.1598230199864138], 'sigma': 0.00726223517112241, 'potential_mean': -2.685551229084199, 'missing_rate': -7.3200624177720595e-06, 'density_quadrature_order': 16}

class PotentialDensityTests(unittest.TestCase):
    def test_saved_cell_retains_original_moments(self):
        f=FIXTURE;r=np.asarray(f['raw']);poly=np.asarray(f['polygon']);m=r[0]
        target=m*f['potential_mean']+.00025*f['missing_rate']
        with patch.dict('os.environ',LAGFORCE_DENSITY_QUADRATURE_ORDER='16'):
            d=PotentialDensity(poly,m,r[1:3],r[3:7].reshape(2,2),target,phi,lambda p:polar_rule(p,16,radii=(.25,.5)))
        x,w=polar_rule(poly,32,radii=(.25,.5));dm=w*d.value(x)
        self.assertLess(abs(dm.sum()/m-1),1e-11)
        np.testing.assert_allclose(dm@x,r[1:3],atol=1e-12,rtol=1e-11)
        np.testing.assert_allclose(x.T@(dm[:,None]*x),r[3:7].reshape(2,2),atol=1e-12,rtol=1e-11)
        self.assertLess(abs(dm@phi(x)-target),1e-12)
        self.assertTrue(np.all(d.value(poly)>0))
        np.testing.assert_allclose(d.weights*d.value(d.points),d.weighted_mass,atol=1e-17,rtol=1e-12)

    def test_dependent_potential_rejected(self):
        f=FIXTURE;r=np.asarray(f['raw']);poly=np.asarray(f['polygon'])
        with self.assertRaises(ValueError):
            PotentialDensity(poly,r[0],r[1:3],r[3:7].reshape(2,2),r[0],lambda x:np.ones(len(x)),lambda p:polar_rule(p,16))

if __name__=='__main__':unittest.main()

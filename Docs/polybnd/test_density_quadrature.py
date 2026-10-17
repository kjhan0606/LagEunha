import os, unittest
from unittest.mock import patch
import numpy as np
from moment_ale_evolve import decode
from density_quadrature import rule, validate_order
from periodic2d_hllc import quadrature
from physical_field_accuracy_audit import environment
from kepler_disk_reference import polar_rule

FIXTURE = {'polygon': [[1.0388127391164852, 3.886931494252504], [1.033745768210947, 3.884960100963199], [0.963203699739754, 3.718523896829917], [0.9697817868818536, 3.7037115209079006], [1.1361488987625854, 3.6460366892185823], [1.2061796513326293, 3.811230394499819]], 'raw': [0.005826010545512806, 0.006718124880253451, 0.021777362181944572, 0.007752219100203492, 0.025116363204068348, 0.025116363204068345, 0.08141659341981328, -0.0015322074582688595, -0.0038344827335800875, -0.0017690294759013137, -0.004421390640964348, -0.005732023203117598, -0.014331043639904475, 0.0004046962253931312, 0.001007811603862542, 0.0010078116038625422, 0.0025240796217565315, 0.001465853419309265]}

class DensityQuadratureTests(unittest.TestCase):
    def test_default_rule_exact(self):
        p=np.asarray(FIXTURE['polygon'])
        x,w=quadrature(p);xx,ww=rule(p,6)
        np.testing.assert_array_equal(x,xx);np.testing.assert_array_equal(w,ww)

    def test_config(self):
        self.assertEqual(environment({})['LAGFORCE_DENSITY_QUADRATURE_ORDER'],'6')
        self.assertEqual(environment({'density_quadrature_order':16})['LAGFORCE_DENSITY_QUADRATURE_ORDER'],'16')
        for value in (0,5,65,12.5,True,'bad'):
            with self.assertRaises((ValueError,TypeError)):validate_order(value)

    def test_saved_cell_independent_mass(self):
        p=np.asarray(FIXTURE['polygon']);raw=np.asarray(FIXTURE['raw']);saved=raw.copy()
        x,w=polar_rule(p,24,radii=(.25,.5));errors=[]
        for order in (6,16):
            with patch.dict(os.environ,{'LAGFORCE_DENSITY_QUADRATURE_ORDER':str(order)}):
                d=decode(p,raw)['density']
            self.assertEqual(d.quadrature_order,order)
            errors.append(abs(w@d.value(x)/raw[0]-1))
        self.assertGreater(errors[0],1e-6);self.assertLess(errors[1],1e-11)
        np.testing.assert_array_equal(raw,saved)

if __name__=='__main__':unittest.main()

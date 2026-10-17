import unittest
import numpy as np
from periodic2d_shear import initial,diagnose,P0
from periodic2d_hllc import rhs,GAMMA


class ShearTests(unittest.TestCase):
    def test_initial_variance_and_control(self):
        m,y,b=initial(4,.5,'conservative');mr,yr,br=initial(4,.5,'resolved')
        np.testing.assert_array_equal(m,mr);np.testing.assert_array_equal(y[:,:5],yr[:,:5])
        self.assertAlmostEqual((y[:,5]-yr[:,5]).sum(),b['unresolved_kinetic'],places=13)
        d,_=diagnose(y,m,.5);dr,_=diagnose(yr,mr,.5)
        self.assertAlmostEqual(d['pressure_mean']-P0,b['pressure_bias_mean'],places=13)
        self.assertLess(d['pressure_rms_from_projected_exact'],1e-13)
        self.assertLess(dr['pressure_rms_from_p0'],1e-13)
        self.assertAlmostEqual(d['total_energy'],P0/(GAMMA-1)+.5*.5**2/2,places=11)

    def test_rhs_conservation_and_velocity(self):
        for closure in ('conservative','resolved'):
            m,y,_=initial(4,.5,closure);dy,d=rhs(y,m)
            self.assertLess(d['mass_flux'],2e-13)
            np.testing.assert_allclose(dy[:,3:].sum(0),0.,atol=1e-14)
            np.testing.assert_array_equal(dy[:,:2],y[:,3:5]/m[:,None])


if __name__=='__main__':unittest.main()

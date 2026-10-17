import unittest
import numpy as np
from scipy.integrate import quad
from kepler_spline_gravity import acceleration,potential,enclosed_mass


class SplineGravityTests(unittest.TestCase):
    def test_potential_matches_independent_force_integral(self):
        for r in (0.,.01,.1,.249,.25,.251,.4,.5,.7,2.):
            expected=-2.-quad(lambda s: float(enclosed_mass(np.array(s/.5)))/s**2 if s else 0.,r,.5,points=[.25] if r<.25 else None,epsabs=1e-12)[0] if r<.5 else -1/r
            self.assertAlmostEqual(float(potential([r,0.])),expected,places=12)

    def test_acceleration_is_negative_potential_gradient(self):
        x=np.array([[.03,.04,.02],[.12,-.2,.03],[.26,.04,.1],[.51,.2,-.1],[1.,2.,.3]])
        h=2e-6;g=np.zeros_like(x)
        for j in range(3):
            delta=np.eye(3)[j]*h;g[:,j]=-(potential(x+delta)-potential(x-delta))/(2*h)
        np.testing.assert_allclose(g,acceleration(x),atol=3e-8,rtol=0)

    def test_origin_and_newtonian_exterior(self):
        np.testing.assert_array_equal(acceleration([0.,0.]),[0.,0.])
        self.assertAlmostEqual(float(potential([0.,0.])), -5.6,places=13)
        x=np.array([[1.,0.],[0.,2.],[3.,4.]])
        np.testing.assert_allclose(acceleration(x),-x/np.linalg.norm(x,axis=1)[:,None]**3,atol=1e-15,rtol=0)

    def test_center_translation_and_invalid_parameters(self):
        x=np.array([[.2,.3],[.4,-.1]]);c=np.array([2.,-3.])
        np.testing.assert_allclose(acceleration(x+c,center=c),acceleration(x),atol=2e-14,rtol=0)
        with self.assertRaises(ValueError):potential(x,rc=0.)
        with self.assertRaises(ValueError):acceleration(x,center=[0.,0.,0.])

if __name__=='__main__':unittest.main()

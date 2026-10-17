import unittest
import numpy as np
from positive_exponential_tilt import fit


class PositiveTiltTests(unittest.TestCase):
    def test_known_positive_tilt(self):
        x,w=np.polynomial.legendre.leggauss(24);F=np.column_stack([x,x*x,x**3])
        known=np.array([.4,-.2,.1]);p=w*np.exp(F@known);p/=p.sum();target=p@F;old=target.copy()
        theta,logZ,actual,diag=fit(w,F,target)
        np.testing.assert_allclose(theta,known,atol=2e-10,rtol=0)
        np.testing.assert_allclose(actual@F,target,atol=2e-12,rtol=0)
        np.testing.assert_array_equal(target,old);self.assertTrue(np.all(actual>0))

    def test_invalid_weights_rejected(self):
        with self.assertRaises(ValueError):fit([1.,-1.],[[0.],[1.]],[.5])

if __name__=='__main__':unittest.main()

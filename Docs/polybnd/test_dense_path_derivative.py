import unittest
import numpy as np
from numpy.polynomial import Polynomial as P
from scipy.integrate._ivp.rk import Dop853DenseOutput
from dense_path_derivative import segment_derivative
class DenseTests(unittest.TestCase):
    def test_exact_polynomial_derivative(self):
        F=np.random.default_rng(417).normal(size=(7,3));old=np.array([.2,-.3,.8]);seg=Dop853DenseOutput(.2,.9,old,F)
        polynomials=[]
        for k in range(3):
            p=P([0.])
            for i,f in enumerate(reversed(F)):p=(p+f[k])*(P([0,1]) if i%2==0 else P([1,-1]))
            polynomials.append(p+old[k])
        for t in (.2,.35,.7,.9):
            x=(t-.2)/.7
            np.testing.assert_allclose(seg(t),[p(x) for p in polynomials],atol=2e-14,rtol=0)
            np.testing.assert_allclose(segment_derivative(seg,t),[p.deriv()(x)/.7 for p in polynomials],atol=3e-14,rtol=0)
if __name__=='__main__':unittest.main()

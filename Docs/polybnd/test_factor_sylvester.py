import unittest
import numpy as np
from factor_moment_evolve import advance_factors_sylvester
class SylvesterTests(unittest.TestCase):
    def test_loss_tangent_full_and_deficient_rank(self):
        K=np.array([[2.,.4],[.4,1.]])
        for L in [np.array([[.3,0.],[.1,.2]]),np.array([[.3,.1],[.6,.2]]),np.zeros((2,2))]:
            d={'cov_A':[np.zeros((2,2))],'cov_K':[K],'cov_Qfactor':[np.zeros((2,0))]}
            dt=1e-7;M=advance_factors_sylvester([L],d,dt)[0]
            np.testing.assert_allclose((M@M.T-L@L.T)/dt,-L@K@L.T,atol=3e-7,rtol=0)
    def test_constant_source_loss_local_order_three(self):
        errors=[];k=2.;q=.3;r0=.2
        d={'cov_A':[np.zeros((2,2))],'cov_K':[k*np.eye(2)],'cov_Qfactor':[np.sqrt(q)*np.eye(2)]}
        for dt in (.02,.01,.005):
            L=advance_factors_sylvester([np.sqrt(r0)*np.eye(2)],d,dt)[0]
            exact=(np.exp(-k*dt)*r0+q*(-np.expm1(-k*dt))/k)*np.eye(2)
            errors.append(np.linalg.norm(L@L.T-exact))
            self.assertGreater(np.linalg.eigvalsh(L@L.T).min(),0)
        for a,b in zip(errors[:-1],errors[1:]):self.assertGreater(np.log2(a/b),2.95)
if __name__=='__main__':unittest.main()

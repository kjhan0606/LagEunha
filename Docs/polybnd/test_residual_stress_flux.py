import unittest
import numpy as np
from compatible_moment_rhs import gaussian_flux
class ResidualStressTests(unittest.TestCase):
    def test_positive_velocity_quadrature(self):
        mean=np.array([.3,-.2]);R=np.array([[.04,.01],[.01,.03]]);rho=1.7
        n=np.array([.6,.8]);wn=.07;L=np.linalg.cholesky(R)
        velocities=np.array([mean+s*np.sqrt(2)*L[:,k] for k in (0,1) for s in (-1,1)])
        rates=rho*(velocities@n-wn)/4;fm=sum(rates)
        p=rates@velocities;Pi=velocities.T@(rates[:,None]*velocities)
        E=rates@(.5*np.sum(velocities*velocities,axis=1))
        result=gaussian_flux(mean,R,rho,n,fm)
        for x,y in zip(result,(p,Pi,E)):np.testing.assert_allclose(x,y,atol=2e-16,rtol=0)
if __name__=='__main__':unittest.main()

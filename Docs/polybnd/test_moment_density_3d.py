import unittest
import numpy as np
from moment_density_3d import quadrature,Density3D,features

CUBE=np.array([[i,j,k] for i in (0.,1.) for j in (0.,1.) for k in (0.,1.)])


class Density3DTests(unittest.TestCase):
    def test_positive_quadrature_matches_cube_momials(self):
        x,w=quadrature(CUBE,5)
        self.assertGreater(w.min(),0.)
        for a,b,c in ((0,0,0),(1,0,0),(2,1,0),(2,2,2),(3,1,1)):
            self.assertAlmostEqual(w@(x[:,0]**a*x[:,1]**b*x[:,2]**c),1/((a+1)*(b+1)*(c+1)),places=13)

    def test_uniform_density_and_affine_cell(self):
        A=np.array([[1.2,.3,.2],[.1,.8,-.1],[0.,.2,1.1]]);offset=np.array([2.,-1.,3.]);vertices=CUBE@A.T+offset
        x,w=quadrature(vertices);rho=1.7;m=rho*w.sum();q=rho*(w@x);Q=rho*x.T@(w[:,None]*x)
        d=Density3D(vertices,m,q,Q)
        np.testing.assert_allclose(d.value(vertices),rho,rtol=0,atol=2e-11)
        self.assertAlmostEqual(d.volume,abs(np.linalg.det(A)),places=13)

    def test_nonuniform_moments_with_independent_quadrature(self):
        x,w=quadrature(CUBE,10);theta=np.array([.3,-.2,.1,.4,.1,-.1,.2,.05,-.3]);rho=np.exp(features(x,np.zeros(3),1.)@theta)
        dm=w*rho;m=dm.sum();q=dm@x;Q=x.T@(dm[:,None]*x)
        d=Density3D(CUBE,m,q,Q,order=8)
        fine_x,fine_w=quadrature(CUBE,11);fine_dm=fine_w*d.value(fine_x)
        self.assertGreater(d.value(CUBE).min(),0.)
        self.assertAlmostEqual(fine_dm.sum(),m,places=12)
        np.testing.assert_allclose(fine_dm@fine_x,q,atol=2e-12,rtol=0)
        np.testing.assert_allclose(fine_x.T@(fine_dm[:,None]*fine_x),Q,atol=2e-12,rtol=0)
        np.testing.assert_allclose(d.value(fine_x),np.exp(features(fine_x,np.zeros(3),1.)@theta),atol=2e-10,rtol=0)

    def test_invalid_covariance_and_infeasible_center_rejected(self):
        with self.assertRaisesRegex(ValueError,'covariance'):Density3D(CUBE,1.,np.ones(3),np.eye(3)*.1)
        q=np.array([2.,.5,.5]);Q=np.outer(q,q)+np.eye(3)/12
        with self.assertRaises(RuntimeError):Density3D(CUBE,1.,q,Q)

if __name__=='__main__':unittest.main()

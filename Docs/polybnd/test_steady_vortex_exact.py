import unittest
import numpy as np
from steady_vortex_exact import fields,parameters,RADIUS,CENTER,physical_raw


class SteadyVortexTests(unittest.TestCase):
    def test_stationary_euler_residual_and_positive_pressure(self):
        g=np.linspace(.02,.98,25);x,y=np.meshgrid(g,g);points=np.column_stack([x.ravel(),y.ravel()])
        rho,v,P,B,gradient=fields(points,derivatives=True)
        np.testing.assert_allclose(np.trace(B,axis1=-2,axis2=-1),0.,atol=2e-15)
        np.testing.assert_allclose(np.einsum('nij,nj->ni',B,v)+gradient,0.,atol=5e-15)
        np.testing.assert_allclose(np.einsum('ni,ni->n',v,gradient),0.,atol=2e-15)
        self.assertGreater(P.min(),.6);np.testing.assert_array_equal(rho,np.ones(len(points)))

    def test_analytic_gradient_against_finite_difference(self):
        points=np.array([[.4,.52],[.65,.6],[.1,.1],[.5,.5]]);eps=2e-5
        _,v,P,B,g=fields(points,derivatives=True)
        for j in range(2):
            delta=np.eye(2)[j]*eps;_,vp,pp=fields(points+delta);_,vm,pm=fields(points-delta)
            _,vpp,ppp=fields(points+2*delta);_,vmm,pmm=fields(points-2*delta)
            np.testing.assert_allclose((-vpp+8*vp-8*vm+vmm)/(12*eps),B[:,:,j],atol=5e-10,rtol=0.)
            np.testing.assert_allclose((-ppp+8*pp-8*pm+pmm)/(12*eps),g[:,j],atol=1e-10,rtol=0.)

    def test_peak_and_reflecting_outer_wall(self):
        radius=RADIUS/np.sqrt(7);_,v,_=fields([CENTER+[radius,0.]])
        self.assertAlmostEqual(np.linalg.norm(v),.5,places=14)
        _,v,P=fields([[0.,.5],[1.,.5],[.5,0.],[.5,1.]])
        np.testing.assert_array_equal(v,np.zeros((4,2)));np.testing.assert_array_equal(P,np.ones(4))
        with self.assertRaises(ValueError):parameters(2.)

    def test_positive_raw_quadrature_and_global_symmetry(self):
        # Four patches avoid relying on a single full-box curved-support rule.
        raw=np.zeros(18)
        for i in range(4):
            for j in range(4):
                p=np.array([[i,j],[i+1,j],[i+1,j+1],[i,j+1]])/4
                r,points,w,*_=physical_raw(p,0.,order=24);raw+=r
                self.assertTrue(np.all(w>0))
        self.assertAlmostEqual(raw[0],1.,places=13)
        np.testing.assert_allclose(raw[7:9],0.,atol=2e-15)
        self.assertGreater(raw[10]-raw[11],0.)
        self.assertGreater(raw[-1],1.)


if __name__=='__main__':unittest.main()

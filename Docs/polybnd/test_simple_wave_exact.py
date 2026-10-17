import unittest
import numpy as np
from simple_wave_exact import fields,physical_raw,GAMMA,P0,C0


class SimpleWaveExactTests(unittest.TestCase):
    def test_nonlinear_euler_equations_before_shock(self):
        x=np.linspace(.05,.95,137);points=np.column_stack([x,np.full(len(x),.37)]);t=.07;h=2e-5
        def conserved(p,time):
            rho,v,P=fields(p,time,.02);u=v[:,0];E=P/(GAMMA-1)+rho*u*u/2
            return np.column_stack([rho,rho*u,E]),np.column_stack([rho*u,rho*u*u+P,(E+P)*u])
        time_difference=(conserved(points,t+h)[0]-conserved(points,t-h)[0])/(2*h)
        dx=np.array([h,0.]);space_difference=(conserved(points+dx,t)[1]-conserved(points-dx,t)[1])/(2*h)
        self.assertLess(float(np.max(abs(time_difference+space_difference))),3e-7)
        rho,v,P=fields(points,t,.02)
        np.testing.assert_allclose(P/rho**GAMMA,P0,atol=2e-15)
        c=np.sqrt(GAMMA*P/rho)
        np.testing.assert_allclose(v[:,0]-2*c/(GAMMA-1),-2*C0/(GAMMA-1),atol=2e-15)

    def test_reference_guards_and_fixed_boundary(self):
        points=np.array([[0.,0.],[1.,.4],[0.,1.]])
        rho,v,P=fields(points,.1,.02)
        np.testing.assert_array_equal(v,np.zeros_like(v));np.testing.assert_array_equal(rho,np.ones(3))
        with self.assertRaises(ValueError):fields(points,1.,.02)
        with self.assertRaises(ValueError):fields(points,.2,2.)

    def test_high_order_moments_have_positive_velocity_gram(self):
        poly=np.array([[.23,.2],[.37,.2],[.37,.7],[.23,.7]])
        r,*_=physical_raw(poly,0.,.02,12);fine,*_=physical_raw(poly,0.,.02,16)
        np.testing.assert_allclose(r,fine,rtol=1e-12,atol=1e-14)
        S=r[3:7].reshape(2,2)-np.outer(r[1:3],r[1:3])/r[0]
        C=r[9:13].reshape(2,2)-np.outer(r[1:3],r[7:9])/r[0]
        B=np.linalg.solve(S,C).T
        R=r[13:17].reshape(2,2)-np.outer(r[7:9],r[7:9])/r[0]-B@S@B.T
        self.assertGreaterEqual(np.linalg.eigvalsh(R).min(),-1e-16)

    def test_analytic_time_derivatives_match_centered_difference(self):
        points=np.column_stack([np.linspace(.1,.9,103),np.full(103,.3)]);t=.07;h=1e-6
        values=fields(points,t,.02,True);plus=fields(points,t+h,.02);minus=fields(points,t-h,.02)
        for derivative,a,b in zip(values[3:],plus,minus):
            np.testing.assert_allclose(derivative,(a-b)/(2*h),rtol=1e-5,atol=2e-9)


if __name__=='__main__':unittest.main()

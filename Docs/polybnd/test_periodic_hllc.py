import unittest
import numpy as np
from periodic_hllc_convergence import geometry,initial,rhs,exact,GAMMA


class PeriodicTests(unittest.TestCase):
    def test_geometry_and_face_velocity(self):
        x=np.array([.08,.32,.57,.84]); w=np.array([.001,-.002,.0015,0.])
        v=np.array([.1,-.03,.04,.07]); z=np.array([.003,-.001,0.,.002])
        l,r,a,d=geometry(x,w); self.assertAlmostEqual(a.sum(),1.)
        speed=(v+np.roll(v,-1))/2+(z-np.roll(z,-1))/(2*d)-(w-np.roll(w,-1))*(np.roll(v,-1)-v)/(2*d*d)
        eps=1e-6
        fd=(geometry(x+eps*v,w+eps*z)[1]-geometry(x-eps*v,w-eps*z)[1])/(2*eps)
        np.testing.assert_allclose(fd,speed,atol=1e-10)
        np.testing.assert_allclose(geometry(x,w+1)[2],a,atol=1e-14)

    def test_translating_contact_rhs(self):
        for kind in ('contact','smooth_contact'):
            m,y=initial(16,kind,.7); dy,d=rhs(y,m)
            np.testing.assert_allclose(dy[:,0],.7,atol=1e-14)
            np.testing.assert_allclose(dy[:,1:],0.,atol=2e-13)

    def test_periodic_conservation(self):
        m,y=initial(16,'wave',.7); dy,d=rhs(y,m)
        self.assertLess(d['mass_flux'],2e-13)
        np.testing.assert_allclose(dy[:,2:].sum(0),0.,atol=1e-14)
        np.testing.assert_array_equal(dy[:,0],y[:,2]/m)

    def test_exact_euler_residual(self):
        x=np.linspace(.02,.98,20); t=.2; h=1e-5
        def f(x,t):
            u=exact(x,t,'wave',.7); r,p,e=u.T; v=p/r
            pressure=(GAMMA-1)*(e-.5*p*v)
            return np.column_stack([p,p*v+pressure,(e+pressure)*v])
        ut=(exact(x,t+h,'wave',.7)-exact(x,t-h,'wave',.7))/(2*h)
        fx=(f(x+h,t)-f(x-h,t))/(2*h)
        np.testing.assert_allclose(ut+fx,0.,atol=1e-9)


if __name__=='__main__': unittest.main()

import unittest
import numpy as np
from fixedmass_hllc import fan,sample,state,rhs,initial_case


class HLLCTests(unittest.TestCase):
    def test_contact_moving_and_stationary(self):
        n=np.array([.6,.8])
        for speed in (0.,.7,-.9):
            v=speed*n
            f=fan(state(1.,v,1.),state(100.,v,1.),n)
            flux,_,_=sample(f,speed)
            np.testing.assert_allclose(flux,np.r_[0.,n,speed],atol=2e-13)

    def test_consistent_flux(self):
        n=np.array([1.,0.]); u=state(2.,[.3,.2],.9)
        for wn in (-3.,0.,.2,.3,.7,3.):
            flux,_,_=sample(fan(u,u,n),wn)
            exact=np.r_[u[0]*.3,u[1:3]*.3+.9*n,(u[3]+.9)*.3]-wn*u
            np.testing.assert_allclose(flux,exact,atol=2e-14)

    def test_reverse_and_slope(self):
        n=np.array([.6,.8]); l=state(1.,[.2,-.1],1.); r=state(.7,[-.1,.3],.8)
        f=fan(l,r,n); rev=fan(r,l,-n)
        for wn in (-3.,-.2,.4,3.):
            flux,rho,_=sample(f,wn)
            np.testing.assert_allclose(flux,-sample(rev,-wn)[0],atol=1e-14)
            slope=(sample(f,wn+1e-6)[0][0]-sample(f,wn-1e-6)[0][0])/2e-6
            self.assertAlmostEqual(slope,-rho,places=8)

    def test_static_contact(self):
        for kind in ('uniform','contact','contact100'):
            m,y=initial_case(kind); dy,d=rhs(y,m)
            np.testing.assert_allclose(dy,0.,atol=2e-13)

    def test_conservation_and_no_steering(self):
        for kind in ('acoustic','mixed'):
            m,y=initial_case(kind); dy,d=rhs(y,m)
            self.assertLess(d['mass_flux'],2e-13)
            self.assertLess(abs(dy[20:24].sum()),2e-13)
            np.testing.assert_allclose(dy[12:20].reshape(4,2).sum(0),dy[24:26],atol=2e-13)
            np.testing.assert_array_equal(dy[:8],y[12:20]/np.repeat(m,2))


if __name__=='__main__': unittest.main()

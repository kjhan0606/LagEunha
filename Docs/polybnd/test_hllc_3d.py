import unittest
import numpy as np
from fixedmass_hllc import state,fan,sample,GAMMA


class HLLC3DTests(unittest.TestCase):
    def test_moving_contact_has_zero_material_mass_flux(self):
        n=np.array([1.,2.,3.]);n/=np.linalg.norm(n);v=np.array([.2,-.3,.1]);P=.7
        result=fan(state(1.,v,P),state(.125,v,P),n)
        flux,_,_=sample(result,v@n)
        np.testing.assert_allclose(flux,np.r_[0.,P*n,P*(v@n)],atol=3e-15,rtol=0.)

    def test_three_dimensional_rotation_covariance(self):
        rng=np.random.default_rng(73);R,_=np.linalg.qr(rng.normal(size=(3,3)))
        n=np.array([.3,.4,np.sqrt(.75)]);vl=np.array([.2,.1,-.3]);vr=np.array([-.1,.3,.2])
        old=fan(state(1.,vl,1.),state(.4,vr,.3),n)
        new=fan(state(1.,R@vl,1.),state(.4,R@vr,.3),R@n)
        np.testing.assert_allclose(new[0],old[0],atol=2e-15)
        for k in (1,2):
            expected=old[k].copy();expected[:,1:-1]=expected[:,1:-1]@R.T
            np.testing.assert_allclose(new[k],expected,atol=3e-15)

    def test_ale_galilean_transform(self):
        n=np.array([.6,0.,.8]);vl=np.array([.4,.2,-.1]);vr=np.array([-.2,.1,.3]);boost=np.array([.13,-.21,.17]);wn=.07
        f,_,_=sample(fan(state(1.,vl,1.),state(.3,vr,.2),n),wn)
        g,_,_=sample(fan(state(1.,vl+boost,1.),state(.3,vr+boost,.2),n),wn+boost@n)
        expected=np.r_[f[0],f[1:-1]+boost*f[0],f[-1]+boost@f[1:-1]+.5*(boost@boost)*f[0]]
        np.testing.assert_allclose(g,expected,atol=3e-15)

    def test_uniform_state_and_dimension_guard(self):
        n=np.array([0.,0.,1.]);v=np.array([.2,.3,.4]);u=state(.7,v,.5)
        f,_,_=sample(fan(u,u,n),0.)
        np.testing.assert_allclose(f,np.r_[.7*.4,.7*v*.4+.5*n,(u[-1]+.5)*.4],atol=1e-15)
        with self.assertRaises(ValueError):fan(u,u,n[:2])

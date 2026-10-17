import unittest
import numpy as np
from moment_density_3d import quadrature
from moment_gravity_kick import moments
from moment_state_3d import pack,unpack,decode,face_flux,ALPHA

CUBE=np.array([[i,j,k] for i in (0.,1.) for j in (0.,1.) for k in (0.,1.)])


class State3DTests(unittest.TestCase):
    def test_three_dimensional_factor_injection_and_rank_deficient_loss(self):
        from factor_moment_evolve import advance_factors_sylvester
        dt=.03;L=np.zeros((3,3));F=np.diag([.1,.2,.3])
        d=dict(cov_A=[np.zeros((3,3))],cov_K=[np.eye(3)],cov_Qfactor=[F])
        result=advance_factors_sylvester([L],d,dt)[0]
        np.testing.assert_allclose(result@result.T,dt*F@F.T,atol=1e-16,rtol=0)
        L=np.diag([.1,0.,.3]);K=np.diag([2.,1.,3.])
        d.update(cov_K=[K],cov_Qfactor=[np.zeros((3,0))])
        result=advance_factors_sylvester([L],d,dt)[0]
        expected=np.diag(np.diag(L)**2*np.exp(-dt*np.diag(K)))
        np.testing.assert_allclose(result@result.T,expected,atol=1e-16,rtol=0)

    def test_affine_velocity_and_positive_covariance_recovered(self):
        x,w=quadrature(CUBE);dm=1.2*w;B=np.array([[.1,.2,-.1],[-.3,.05,.2],[.1,0.,.2]]);v0=np.array([.2,-.1,.3]);L=np.diag([.03,.02,.01])
        v=(v0+x@B.T)[:,None,:]+(ALPHA@L.T)[None,:,:]
        raw=moments(x,dm,v,2.);r=pack(*(raw[k] for k in ('m','q','Q','p','D','Pi','E')))
        self.assertEqual(r.shape,(35,));s=decode(CUBE,r)
        np.testing.assert_allclose(s['B'],B,atol=1e-14,rtol=0)
        np.testing.assert_allclose(s['R'],dm.sum()*L@L.T,atol=1e-14,rtol=0)
        self.assertAlmostEqual(s['P'],4/3,places=13)
        np.testing.assert_allclose(s['v'],v0+B@np.full(3,.5),atol=1e-14,rtol=0)

    def test_uniform_face_matches_euler_flux(self):
        x=np.array([.2,.3,.4]);n=np.array([.6,0.,.8]);v=np.array([.3,-.2,.1]);rho=1.4;P=.7;wn=.05;s=v@n-wn
        left=np.tile(v,(6,1));right=left.copy();a=np.full(6,rho*max(s,0)/6);b=np.full(6,rho*min(s,0)/6)
        f=face_flux(x,n,left,right,a,b,P/((5/3-1)*rho),P/((5/3-1)*rho),P,v@n)
        m,q,Q,p,D,Pi,E=unpack(f)
        self.assertAlmostEqual(m,rho*s,places=14)
        np.testing.assert_allclose(p,rho*s*v+P*n,atol=1e-14,rtol=0)
        np.testing.assert_allclose(D,np.outer(x,p),atol=1e-14,rtol=0)
        self.assertAlmostEqual(E,s*(.5*rho*(v@v)+P/(5/3-1))+P*(v@n),places=14)

    def test_face_reversal_has_equal_opposite_all_moments(self):
        rng=np.random.default_rng(129);vl=rng.normal(size=(6,3));vr=rng.normal(size=(6,3));a=rng.uniform(.01,.2,6);b=-rng.uniform(.01,.2,6)
        x=np.array([.2,.3,.4]);n=np.array([0.,0.,1.])
        f=face_flux(x,n,vl,vr,a,b,1.,2.,.7,.2)
        g=face_flux(x,-n,vr,vl,-b,-a,2.,1.,.7,-.2)
        np.testing.assert_allclose(f,-g,atol=2e-15,rtol=0)

if __name__=='__main__':unittest.main()

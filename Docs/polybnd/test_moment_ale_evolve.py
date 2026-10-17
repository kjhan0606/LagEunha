import unittest
import numpy as np
from moment_ale_evolve import Density,features,initial,rhs,geometry,decode,P0
from periodic2d_hllc import quadrature


class MomentEvolutionTests(unittest.TestCase):
    def setUp(self):
        import moment_ale_evolve as model
        from adaptive_moment_batch import install_affine_shear_domain
        self.addCleanup(setattr, model, 'domain', model.domain)
        self.addCleanup(setattr, model, 'boundary_velocity', model.boundary_velocity)
        install_affine_shear_domain()
    def test_positive_density_fit(self):
        poly=np.array([[0.,0.],[.8,.1],[1.,.9],[.1,1.]])
        q,w=quadrature(poly);origin=w@q/w.sum();scale=np.sqrt(w.sum())
        rho=np.exp(features(q,origin,scale)@np.array([.2,-.1,.08,.03,-.04]))
        weights=w*rho;m=weights.sum();first=weights@q;second=(q*weights[:,None]).T@q
        fitted=Density(poly,m,first,second)
        self.assertLess(fitted.residual,3e-13)
        np.testing.assert_allclose(fitted.value(q),rho,rtol=1e-10,atol=1e-12)
        self.assertGreater(fitted.value(poly).min(),0.)

    def test_affine_initial_rhs(self):
        mass,y=initial(3);dy,d=rhs(0.,y,mass);state=y[:-4].reshape(9,20);rate=dy[:-4].reshape(9,20)
        np.testing.assert_allclose(rate[:,0],.2+1.1*state[:,1],atol=1e-12)
        np.testing.assert_allclose(rate[:,1],0.,atol=1e-12)
        self.assertLess(d['mass_flux'],2e-13)
        self.assertLess(d['residual_covariance_max'],1e-12)
        np.testing.assert_allclose(rate[:,9:11].sum(0),dy[-4:-2],atol=1e-12)
        self.assertLess(abs(rate[:,19].sum()-dy[-2]),1e-12)
        polys,_=geometry(state[:,:2],state[:,2],0.)
        for p,m,row in zip(polys,mass,state):self.assertAlmostEqual(decode(p,np.r_[m,row[3:]])['P'],P0,places=12)


if __name__=='__main__':unittest.main()

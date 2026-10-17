import unittest
from unittest.mock import patch
import numpy as np
from adaptive_moment_batch import setup
from simple_wave_evolution import install_domain
from factor_moment_evolve import rebuild,advance_factors_sylvester as advance
from kinetic_factor_rhs import rhs
from geometry_mass_lift import rates
from extended_dense_moments import dense_state
import moment_ale_evolve as model


class ExtendedDenseTests(unittest.TestCase):
    def setUp(self):
        # Domain installers mutate module globals; restore them even on failure.
        import moment_ale_evolve as model
        self.addCleanup(setattr, model, 'domain', model.domain)
        self.addCleanup(setattr, model, 'boundary_velocity', model.boundary_velocity)

    def test_sub_ulp_polynomial_changes_are_retained(self):
        y=np.zeros(24);y[3:]=1.;d0=np.zeros(24);d0[3]=1e-16;dm=d0.copy();L=[np.zeros((2,2))]
        # A valid cell covariance; test arithmetic through the actual helper.
        mass,state=setup('wave',1,0.);state=rebuild(state,mass,L);d0=np.zeros_like(state);d0[3]=1e-17
        diag=dict(cov_A=[np.zeros((2,2))],cov_Qfactor=[np.zeros((2,0))])
        path,_,_=dense_state(state,L,d0,d0,diag,mass,1.,.5)
        self.assertEqual(path.dtype,np.dtype(np.longdouble))
        self.assertNotEqual(path[3],np.longdouble(state[3]))

    def test_extended_decode_and_root_preserve_physical_values(self):
        install_domain();mass,state=setup('wave',3,.01);L=[np.zeros((2,2)) for _ in mass];state=rebuild(state,mass,L)
        with patch.dict('os.environ',{'LAGFORCE_VECTOR_KINETIC_ROOT':'1','LAGFORCE_QUADRATIC_VELOCITY':'0'}):
            dy,diag=rhs(0.,state,mass,covariance_factors=L)
            g,gd=rates(0.,state.astype(np.longdouble),L,mass)
        np.testing.assert_allclose(g,dy[:-4].reshape(-1,20)[:,:3],atol=2e-11,rtol=0.)
        self.assertLess(gd['mass_flux'],2e-13)

    def test_extended_density_fit_has_a_smaller_internal_stopping_error(self):
        install_domain();mass,state=setup('wave',3,.01);a=state.astype(np.longdouble)[:-4].reshape(-1,20)
        cells,_=model.geometry(np.asarray(a[:,:2],float),np.asarray(a[:,2],float),0.)
        states=[model.decode(poly,np.r_[m,row[3:]]) for m,row,poly in zip(mass,a,cells)]
        self.assertLessEqual(max(s['density'].residual for s in states),1e-15)

    def test_same_moment_polynomial_and_positive_factors(self):
        install_domain();mass,y=setup('wave',3,.01);L=[np.zeros((2,2)) for _ in mass];y=rebuild(y,mass,L);dt=.0002
        with patch.dict('os.environ',{'LAGFORCE_VECTOR_KINETIC_ROOT':'1','LAGFORCE_QUADRATIC_VELOCITY':'0'}):
            d0,diag0=rhs(0.,y,mass,covariance_factors=L);halfL=advance(L,diag0,dt/2)
            dm,diag=rhs(dt/2,rebuild(y+dt*d0/2,mass,halfL),mass,covariance_factors=halfL)
        from coupled_factor_path_audit import dense_state as original
        before=y.copy()
        for s in [0.,.2,.5,1.]:
            new,factors,_=dense_state(y,L,d0,dm,diag,mass,dt,s);old,oldL,_=original(y,L,d0,dm,diag,mass,dt,s)
            np.testing.assert_allclose(np.asarray(new,float),old,atol=2e-15,rtol=0.)
            for f,g in zip(factors,oldL):np.testing.assert_array_equal(f,g)
        np.testing.assert_array_equal(y,before)


if __name__=='__main__':unittest.main()

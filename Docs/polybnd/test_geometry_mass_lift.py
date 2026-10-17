import unittest
from unittest.mock import patch
import numpy as np
from scipy.integrate import solve_ivp
from geometry_mass_lift import evaluate_polynomial,polynomial_derivative,rates,integrate,lifted_state
from adaptive_moment_batch import setup
from factor_moment_evolve import rebuild,advance_factors_sylvester as advance
from simple_wave_evolution import install_domain,choose_next_dt
from kinetic_factor_rhs import rhs
from smooth_mass_path import evaluate_smooth_path


class GeometryMassLiftTests(unittest.TestCase):
    def setUp(self):
        # Domain installers mutate module globals; restore them even on failure.
        import moment_ale_evolve as model
        self.addCleanup(setattr, model, 'domain', model.domain)
        self.addCleanup(setattr, model, 'boundary_velocity', model.boundary_velocity)

    def test_quadrature_refinement_does_not_limit_physical_step_growth(self):
        checked=dict(ratio=.9,physical_ratio=.05)
        cfg=dict(geometry_step_controller='physical_path')
        self.assertAlmostEqual(choose_next_dt(cfg,.001,.01,checked),.0015)
        self.assertAlmostEqual(choose_next_dt({},.001,.01,checked),.001*.9/np.sqrt(.9))
        self.assertEqual(choose_next_dt(dict(cfg,fixed_dt=True),.001,.01,checked),.01)
        self.assertEqual(choose_next_dt(cfg,.001,.001,checked),.001)

    def test_saved_polynomial_and_its_derivative(self):
        solution=solve_ivp(lambda t,y:np.array([y[1],-y[0]]),(0.,.3),[.2,.4],method='DOP853',rtol=1e-12,atol=1e-13,dense_output=True)
        for v in solution.sol.interpolants:
            segment=dict(t=v.t,t_old=v.t_old,h=v.h,y_old=v.y_old,F=v.F)
            for t in np.linspace(v.t_old,v.t,7):
                np.testing.assert_allclose(evaluate_polynomial(segment,t),solution.sol(t),atol=1e-15,rtol=0.)
                independent=np.imag(evaluate_polynomial(segment,t+1e-30j))/1e-30
                np.testing.assert_allclose(polynomial_derivative(segment,t),independent,atol=2e-14,rtol=0.)

    def test_mass_root_matches_full_rhs_and_independent_audit(self):
        install_domain();mass,y=setup('wave',3,.01);L=[np.zeros((2,2)) for _ in mass];y=rebuild(y,mass,L)
        with patch.dict('os.environ',{'LAGFORCE_VECTOR_KINETIC_ROOT':'1','LAGFORCE_QUADRATIC_VELOCITY':'0'}):
            d,diag=rhs(0.,y,mass,covariance_factors=L);g,gd=rates(0.,y,L,mass)
        np.testing.assert_allclose(g,d[:-4].reshape(-1,20)[:,:3],atol=2e-12,rtol=0.)
        f,audit=evaluate_smooth_path(0.,y,L,mass,(g[:,:2],g[:,2]))
        self.assertLess(np.max(abs(f)),2e-13);self.assertLess(audit['generator_defect'],1e-14)

    def test_geometry_root_can_be_solved_below_the_ode_error_scale(self):
        install_domain();mass,y=setup('wave',3,.01);L=[np.zeros((2,2)) for _ in mass];y=rebuild(y,mass,L)
        with patch.dict('os.environ',{'LAGFORCE_VECTOR_KINETIC_ROOT':'1','LAGFORCE_QUADRATIC_VELOCITY':'0'}):
            old,_=rates(0.,y,L,mass);new,diag=rates(0.,y,L,mass,root_tolerance=1e-16)
        self.assertLess(diag['mass_flux'],1e-16)
        np.testing.assert_allclose(new,old,atol=2e-11,rtol=0.)

    def test_quadratic_mass_root_matches_full_rhs(self):
        from test_quadratic_velocity_reconstruction import fixture
        mass,state,L,*_=fixture()
        with patch.dict('os.environ',{'LAGFORCE_VECTOR_KINETIC_ROOT':'1','LAGFORCE_QUADRATIC_VELOCITY':'1'}):
            d,diag=rhs(0.,state,mass,covariance_factors=L);g,gd=rates(0.,state,L,mass)
            np.testing.assert_allclose(g,d[:-4].reshape(-1,20)[:,:3],atol=2e-12,rtol=0.)
            f,audit=evaluate_smooth_path(0.,state,L,mass,(g[:,:2],g[:,2]))
        self.assertLess(np.max(abs(f)),2e-13);self.assertLess(audit['generator_defect'],1e-14)

    def test_lift_changes_only_geometry_and_keeps_moments_and_factors(self):
        install_domain();mass,y=setup('wave',3,.01);L=[np.zeros((2,2)) for _ in mass];y=rebuild(y,mass,L);dt=.0005
        with patch.dict('os.environ',{'LAGFORCE_VECTOR_KINETIC_ROOT':'1','LAGFORCE_QUADRATIC_VELOCITY':'0'}):
            d0,diag0=rhs(0.,y,mass,covariance_factors=L);halfL=advance(L,diag0,dt/2)
            dm,diag=rhs(dt/2,rebuild(y+dt*d0/2,mass,halfL),mass,covariance_factors=halfL)
            record=dict(t=0.,dt=dt,y=y,L=L,d0=d0,dm=dm,diag={k:diag[k] for k in ['cov_linearization_factors','cov_A','cov_K','cov_Qfactor']})
            saved=y.copy();lift=integrate(record,mass);record['geometry_lift']=lift
            from coupled_factor_path_audit import dense_state
            for s in [0.,.2,.5,.8,1.]:
                path,factors,derivative=lifted_state(record,mass,s)
                original,oldL,_=dense_state(y,L,d0,dm,record['diag'],mass,dt,s)
                np.testing.assert_array_equal(path[:-4].reshape(-1,20)[:,3:],original[:-4].reshape(-1,20)[:,3:])
                np.testing.assert_array_equal(path[-4:],original[-4:])
                for a,b in zip(factors,oldL):np.testing.assert_array_equal(a,b)
                a=derivative[:-4].reshape(-1,20);f,audit=evaluate_smooth_path(s*dt,path,factors,mass,(a[:,:2],a[:,2]))
                self.assertLess(np.max(abs(f)),2e-10);self.assertLess(audit['generator_defect'],2e-8)
            np.testing.assert_array_equal(y,saved)


if __name__=='__main__':unittest.main()

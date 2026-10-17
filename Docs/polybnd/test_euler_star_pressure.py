import unittest
from unittest.mock import patch
import numpy as np
from euler_star_pressure import star


class ExactStarTests(unittest.TestCase):
    def test_sod_agrees_with_separate_reference(self):
        from sod_exact import star as reference
        p,u=star(1.,[0.,0.],1.,.125,[0.,0.],.1,[1.,0.])
        self.assertAlmostEqual(p,reference()['p'],places=13)
        self.assertAlmostEqual(u,reference()['u'],places=13)

    def test_symmetric_rarefaction_has_analytic_positive_pressure(self):
        rho=.025746008938956854;P=.0011513290710691015;v=.17541368706974036;g=5/3;c=np.sqrt(g*P/rho)
        p,u=star(rho,[-v,0.],P,rho,[v,0.],P,[1.,0.])
        expected=P*(1-(g-1)*v/(2*c))**(2*g/(g-1))
        self.assertAlmostEqual(p,expected,places=15);self.assertAlmostEqual(u,0.,places=15)

    def test_normal_reversal_and_galilean_shift_in_3d(self):
        n=np.array([.6,0.,.8]);vl=np.array([.1,.3,.2]);vr=np.array([-.2,.1,.4]);boost=np.array([.3,-.1,.2])
        p,u=star(1.,vl,1.,.2,vr,.3,n);q,w=star(.2,vr,.3,1.,vl,1.,-n)
        self.assertAlmostEqual(p,q,places=14);self.assertAlmostEqual(u,-w,places=14)
        q,w=star(1.,vl+boost,1.,.2,vr+boost,.3,n)
        self.assertAlmostEqual(p,q,places=14);self.assertAlmostEqual(w,u+boost@n,places=14)

    def test_vacuum_is_explicitly_unsupported(self):
        with self.assertRaisesRegex(ValueError,'vacuum'):star(1.,[-10.,0.],1.,1.,[10.,0.],1.,[1.,0.])

    def test_pressure_solver_requires_explicit_driver_configuration(self):
        from simple_wave_evolution import run_with_executor
        with patch.dict('os.environ',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star'),patch('simple_wave_evolution.Path.read_bytes',return_value=b'{}'):
            with self.assertRaisesRegex(ValueError,'pressure Riemann solver'):run_with_executor()

    def test_exact_traction_preserves_covariance_moment_identity(self):
        from test_kinetic_factors import KineticFactorTests
        with patch.dict('os.environ',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star'):
            KineticFactorTests().test_covariance_equation_matches_raw_moments()

    def test_exact_traction_conserves_shared_face_ledgers(self):
        from adaptive_moment_batch import setup
        from factor_moment_evolve import rebuild
        from kinetic_factor_rhs import rhs
        with patch.dict('os.environ',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star'):
            mass,y=setup('wave',3,.01)
            factors=[np.sqrt(m)*np.array([[.01,0.],[.003,.007]]) for m in mass]
            y=rebuild(y,mass,factors)
            dy,diag=rhs(0.,y,mass,covariance_factors=factors)
        rates=dy[:-4].reshape(-1,20)
        totals=np.r_[rates[:,9:11].sum(axis=0),rates[:,19].sum(),
                     (rates[:,12]-rates[:,13]).sum()]
        np.testing.assert_allclose(totals,dy[-4:],atol=1e-13,rtol=0)
        self.assertLess(diag['mass_flux'],2e-13)
        self.assertEqual(diag['generator_velocity_defect'],0.)

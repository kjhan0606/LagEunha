import json,unittest
from unittest.mock import patch
import numpy as np


class VectorMomentFluxTests(unittest.TestCase):
    def test_combined_geometry_horner_and_flux_backends_match_physical_rates(self):
        from test_quadratic_velocity_reconstruction import fixture
        from kinetic_factor_rhs import rhs
        from profile_combined_numerical_kernels import flags,covariance_rates
        mass,state,L,*_=fixture();env=dict(LAGFORCE_ACTIVE_POWER_CELLS='1',LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VECTOR_QUADRATIC_BUDGET='1',LAGFORCE_QUADRATIC_PRESSURE='1',LAGFORCE_PRESSURE_DEGREE='4',LAGFORCE_VECTOR_KINETIC_ROOT='1',LAGFORCE_KINETIC_SMOOTHING='1')
        with patch.dict('os.environ',env),patch.dict('os.environ',flags(False)):old,do=rhs(0.,state,mass,covariance_factors=L)
        with patch.dict('os.environ',env),patch.dict('os.environ',flags(True)):new,dn=rhs(0.,state,mass,covariance_factors=L)
        np.testing.assert_allclose(new,old,atol=3e-12,rtol=1e-10);np.testing.assert_allclose(covariance_rates(do,L),covariance_rates(dn,L),atol=3e-12,rtol=1e-10)

    def test_batched_flux_retains_the_positive_covariance_factor_equation(self):
        from test_quadratic_velocity_reconstruction import QuadraticVelocityTests
        with patch.dict('os.environ',LAGFORCE_VECTOR_MOMENT_FLUX='1'):
            QuadraticVelocityTests().test_positive_factor_equation_matches_raw_covariance_derivative()

    def test_all_primary_and_factor_moment_rates_match_scalar_fluxes(self):
        from test_quadratic_velocity_reconstruction import fixture
        from kinetic_factor_rhs import rhs
        mass,state,L,*_=fixture()
        for quadratic,pressure in [('0','0'),('1','0'),('1','1')]:
            env=dict(LAGFORCE_QUADRATIC_VELOCITY=quadratic,LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VECTOR_QUADRATIC_BUDGET='1',LAGFORCE_QUADRATIC_PRESSURE=pressure,LAGFORCE_PRESSURE_DEGREE='4',LAGFORCE_VECTOR_KINETIC_ROOT='1',LAGFORCE_KINETIC_SMOOTHING='1')
            with patch.dict('os.environ',env),patch.dict('os.environ',LAGFORCE_VECTOR_MOMENT_FLUX='0'):old,do=rhs(0.,state,mass,covariance_factors=L)
            with patch.dict('os.environ',env),patch.dict('os.environ',LAGFORCE_VECTOR_MOMENT_FLUX='1'):new,dn=rhs(0.,state,mass,covariance_factors=L)
            np.testing.assert_allclose(new,old,atol=3e-13,rtol=2e-10)
            for key in ['cov_A','cov_K','cov_Qfactor']:
                for a,b in zip(do[key],dn[key]):np.testing.assert_allclose(a,b,atol=2e-13,rtol=2e-10)
            self.assertEqual(dn['mass_flux'],do['mass_flux']);np.testing.assert_array_equal(new[:-4].reshape(-1,20)[:,:3],old[:-4].reshape(-1,20)[:,:3])

    def test_shared_flux_preserves_global_conservation_ledgers(self):
        from test_quadratic_velocity_reconstruction import fixture
        from kinetic_factor_rhs import rhs
        mass,state,L,*_=fixture();env=dict(LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VECTOR_QUADRATIC_BUDGET='1',LAGFORCE_QUADRATIC_PRESSURE='1',LAGFORCE_PRESSURE_DEGREE='4',LAGFORCE_VECTOR_KINETIC_ROOT='1',LAGFORCE_KINETIC_SMOOTHING='1',LAGFORCE_VECTOR_MOMENT_FLUX='1')
        with patch.dict('os.environ',env):dy,diag=rhs(0.,state,mass,covariance_factors=L)
        a=dy[:-4].reshape(-1,20);np.testing.assert_allclose(a[:,9:11].sum(0),dy[-4:-2],atol=3e-13,rtol=0.);self.assertLess(abs(a[:,19].sum()-dy[-2]),3e-13);self.assertLess(abs(np.sum(a[:,12]-a[:,13])-dy[-1]),3e-13)

    def test_backend_requires_matching_configuration(self):
        from simple_wave_evolution import run_with_executor
        env=dict(LAGFORCE_VECTOR_KINETIC_ROOT='0',LAGFORCE_QUADRATIC_VELOCITY='0',LAGFORCE_VECTOR_QUADRATIC_BUDGET='0',LAGFORCE_QUADRATIC_PRESSURE='0',LAGFORCE_VECTOR_SMOOTH_PATH='0',LAGFORCE_VECTOR_MOMENT_FLUX='0')
        with patch.dict('os.environ',env),patch('simple_wave_evolution.Path.read_bytes',return_value=json.dumps(dict(moment_flux_evaluation_backend='batched')).encode()):
            with self.assertRaisesRegex(ValueError,'moment flux evaluation backend'):run_with_executor()

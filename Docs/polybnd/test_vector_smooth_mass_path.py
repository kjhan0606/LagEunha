import unittest
import json
from unittest.mock import patch
import numpy as np
from smooth_mass_path import evaluate_smooth_path


class VectorSmoothMassPathTests(unittest.TestCase):
    def setUp(self):
        # Domain installers mutate module globals; restore them even on failure.
        import moment_ale_evolve as model
        self.addCleanup(setattr, model, 'domain', model.domain)
        self.addCleanup(setattr, model, 'boundary_velocity', model.boundary_velocity)

    def test_path_backend_configuration_is_explicit(self):
        from simple_wave_evolution import run_with_executor
        cfg=dict(mass_path_evaluation_backend='cell_batched')
        env=dict(LAGFORCE_VECTOR_KINETIC_ROOT='0',LAGFORCE_QUADRATIC_VELOCITY='0',LAGFORCE_VECTOR_QUADRATIC_BUDGET='0',LAGFORCE_QUADRATIC_PRESSURE='0',LAGFORCE_VECTOR_SMOOTH_PATH='0')
        with patch.dict('os.environ',env),patch('simple_wave_evolution.Path.read_bytes',return_value=json.dumps(cfg).encode()):
            with self.assertRaisesRegex(ValueError,'mass path evaluation backend'):run_with_executor()

    def test_all_existing_continuous_path_modes_match_scalar_evaluation(self):
        from test_quadratic_velocity_reconstruction import fixture
        mass,state,L,*_=fixture();rates=np.random.default_rng(443).normal(size=(len(mass),3))*.01
        for quadratic in ['0','1']:
            for coefficient in ['0','1']:
                env=dict(LAGFORCE_QUADRATIC_VELOCITY=quadratic,LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VECTOR_QUADRATIC_BUDGET='1',LAGFORCE_KINETIC_SMOOTHING=coefficient)
                with patch.dict('os.environ',env),patch.dict('os.environ',LAGFORCE_VECTOR_SMOOTH_PATH='0'):old,do=evaluate_smooth_path(0.,state,L,mass,(rates[:,:2],rates[:,2]))
                with patch.dict('os.environ',env),patch.dict('os.environ',LAGFORCE_VECTOR_SMOOTH_PATH='1'):new,dn=evaluate_smooth_path(0.,state,L,mass,(rates[:,:2],rates[:,2]))
                np.testing.assert_allclose(new,old,atol=2e-14,rtol=2e-11);self.assertEqual(dn,do)

    def test_extended_dense_geometry_rates_keep_the_same_arithmetic(self):
        from test_quadratic_velocity_reconstruction import fixture
        mass,state,L,*_=fixture();rates=np.asarray(np.random.default_rng(19).normal(size=(len(mass),3))*.02,dtype=np.longdouble)
        rates+=np.longdouble('1.234567890123456e-18')
        env=dict(LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VECTOR_QUADRATIC_BUDGET='1',LAGFORCE_KINETIC_SMOOTHING='1')
        with patch.dict('os.environ',env),patch.dict('os.environ',LAGFORCE_VECTOR_SMOOTH_PATH='0'):old,do=evaluate_smooth_path(0.,state,L,mass,(rates[:,:2],rates[:,2]))
        with patch.dict('os.environ',env),patch.dict('os.environ',LAGFORCE_VECTOR_SMOOTH_PATH='1'):new,dn=evaluate_smooth_path(0.,state,L,mass,(rates[:,:2],rates[:,2]))
        np.testing.assert_allclose(new,old,atol=2e-14,rtol=2e-11);self.assertEqual(dn,do)

    def test_vector_path_remains_independent_of_instantaneous_mass_root(self):
        from simple_wave_evolution import initial,install_domain
        from kinetic_factor_rhs import rhs
        install_domain();mass,state,L,_=initial(4,.02)
        env=dict(LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VECTOR_QUADRATIC_BUDGET='1',LAGFORCE_KINETIC_SMOOTHING='1',LAGFORCE_VECTOR_KINETIC_ROOT='1',LAGFORCE_VECTOR_SMOOTH_PATH='1')
        with patch.dict('os.environ',env):
            dy,diag=rhs(0.,state,mass,covariance_factors=L);a=dy[:-4].reshape(-1,20)
            with patch('kinetic_factor_rhs.solve',side_effect=RuntimeError('audit must not solve')):flux,detail=evaluate_smooth_path(0.,state,L,mass,(a[:,:2],a[:,2]))
        self.assertLess(float(np.max(abs(flux))),3e-13);self.assertEqual(detail['generator_defect'],0.)

import unittest
from unittest.mock import patch
import os
import numpy as np
from adaptive_moment_batch import setup
from analyze_simple_wave_time import observed,sampled_fields,check_environment


class SimpleWaveTimeAnalysisTests(unittest.TestCase):
    def test_order_omitted_at_the_roundoff_screen(self):
        self.assertAlmostEqual(observed([4e-7,1e-7],1e-12),2.)
        self.assertIsNone(observed([4e-13,1e-13],1e-12))

    def test_common_eulerian_samples_preserve_uniform_fields(self):
        mass,state=setup('wave',3,0.)
        points=np.array([[.12,.13],[.31,.46],[.7,.82],[.92,.2]])
        values=sampled_fields(dict(t=0.,mass=mass,state=state),points)
        np.testing.assert_allclose(values['density'],1.,atol=2e-12)
        np.testing.assert_allclose(values['pressure'],.6,atol=2e-12)
        np.testing.assert_array_equal(values['velocity'],np.zeros_like(points))

    def test_analysis_uses_configured_pressure_degree_and_velocity(self):
        cfg=dict(velocity_reconstruction='quadratic_budget_projected',pressure_reconstruction='polynomial_cell_mean_degree4',quadratic_budget_backend='batched')
        env=dict(LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VECTOR_QUADRATIC_BUDGET='1',LAGFORCE_QUADRATIC_PRESSURE='1',LAGFORCE_PRESSURE_DEGREE='4',LAGFORCE_FAST_PRESSURE_POLYNOMIAL='0')
        with patch.dict(os.environ,env):
            check_environment(cfg)
            with patch.dict(os.environ,LAGFORCE_PRESSURE_DEGREE='2'):
                with self.assertRaisesRegex(ValueError,'pressure_reconstruction'):check_environment(cfg)
            with patch.dict(os.environ,LAGFORCE_QUADRATIC_VELOCITY='0'):
                with self.assertRaisesRegex(ValueError,'velocity_reconstruction'):check_environment(cfg)

    def test_analysis_requires_matching_evaluation_backends(self):
        cfg=dict(velocity_reconstruction='quadratic_budget_projected',pressure_reconstruction='polynomial_cell_mean_degree4',quadratic_budget_backend='batched',pressure_evaluation_backend='horner')
        env=dict(LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VECTOR_QUADRATIC_BUDGET='1',LAGFORCE_QUADRATIC_PRESSURE='1',LAGFORCE_PRESSURE_DEGREE='4',LAGFORCE_FAST_PRESSURE_POLYNOMIAL='1')
        with patch.dict(os.environ,env):
            check_environment(cfg)
            with patch.dict(os.environ,LAGFORCE_FAST_PRESSURE_POLYNOMIAL='0'):
                with self.assertRaisesRegex(ValueError,'pressure_evaluation_backend'):check_environment(cfg)
            with patch.dict(os.environ,LAGFORCE_VECTOR_QUADRATIC_BUDGET='0'):
                with self.assertRaisesRegex(ValueError,'quadratic_budget_backend'):check_environment(cfg)


if __name__=='__main__':unittest.main()

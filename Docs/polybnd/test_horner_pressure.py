import json,unittest
from unittest.mock import patch
import numpy as np
from polynomial_pressure_reconstruction import value,horner_coefficients,prepare


class HornerPressureTests(unittest.TestCase):
    def test_scalar_and_vector_values_are_the_same_polynomial(self):
        rng=np.random.default_rng(11);points=rng.uniform(-1,1,(500,2))
        for degree in [3,4]:
            from polynomial_pressure_reconstruction import powers
            exponent=powers(degree)[1:];meta=dict(degree=degree,mean=.8,center=np.array([.1,-.2]),scale=.7,powers=exponent,basis_means=rng.normal(size=len(exponent)),coefficient=.01*rng.normal(size=len(exponent)))
            reference=value(meta,points);fast=dict(meta,horner_coefficients=horner_coefficients(meta))
            np.testing.assert_allclose(value(fast,points),reference,atol=3e-15,rtol=2e-15)
            for p,r in zip(points[:10],reference[:10]):self.assertAlmostEqual(float(value(fast,p)[0]),r,places=14)

    def test_positivity_decision_and_coefficients_are_identical(self):
        from test_polynomial_pressure_reconstruction import data
        _,_,_,cells,states=data()
        for i,s in enumerate(states):s['P']=1e-5 if i==len(states)//2 else .5+(i%3)
        with patch.dict('os.environ',{'LAGFORCE_FAST_PRESSURE_POLYNOMIAL':'0'}):reference,_=prepare(states,cells,4)
        with patch.dict('os.environ',{'LAGFORCE_FAST_PRESSURE_POLYNOMIAL':'1'}):fast,_=prepare(states,cells,4)
        for a,b,p in zip(reference,fast,cells):
            np.testing.assert_array_equal(a['coefficient'],b['coefficient']);self.assertEqual(a['positivity_scale'],b['positivity_scale']);np.testing.assert_allclose(value(a,p),value(b,p),atol=3e-13,rtol=2e-13)

    def test_runtime_backend_must_match_the_configuration(self):
        from simple_wave_evolution import run_with_executor
        cfg=dict(mass_root_backend='legacy_dense',velocity_reconstruction='affine',pressure_reconstruction='polynomial_cell_mean_degree4',pressure_evaluation_backend='horner')
        flags=dict(LAGFORCE_QUADRATIC_VELOCITY='0',LAGFORCE_QUADRATIC_PRESSURE='1',LAGFORCE_PRESSURE_DEGREE='4',LAGFORCE_VECTOR_QUADRATIC_BUDGET='0',LAGFORCE_VECTOR_KINETIC_ROOT='0',LAGFORCE_FAST_PRESSURE_POLYNOMIAL='0')
        with patch.dict('os.environ',flags),patch('simple_wave_evolution.Path.read_bytes',return_value=json.dumps(cfg).encode()):
            with self.assertRaisesRegex(ValueError,'pressure evaluation backend'):run_with_executor()

    def test_time_analysis_rejects_a_different_field_model(self):
        from analyze_simple_wave_time import check_environment
        flags=dict(LAGFORCE_QUADRATIC_VELOCITY='0',LAGFORCE_QUADRATIC_PRESSURE='1',LAGFORCE_PRESSURE_DEGREE='4',LAGFORCE_VECTOR_QUADRATIC_BUDGET='0',LAGFORCE_FAST_PRESSURE_POLYNOMIAL='0')
        with patch.dict('os.environ',flags):
            check_environment(dict(pressure_reconstruction='polynomial_cell_mean_degree4'))
            with self.assertRaisesRegex(ValueError,'environment mismatch'):check_environment(dict(pressure_reconstruction='linear'))


if __name__=='__main__':unittest.main()

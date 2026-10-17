import unittest
import os
from unittest.mock import patch
import numpy as np
from compare_physical_budget import scales_and_goals,compatible,field_difference,endpoint_matches,sampled_pair


class PhysicalBudgetComparisonTests(unittest.TestCase):
    def test_endpoint_comparison_allows_only_time_roundoff(self):
        self.assertTrue(endpoint_matches(np.nextafter(.05,0.),.05));self.assertFalse(endpoint_matches(.05-1e-8,.05));self.assertFalse(endpoint_matches(float('nan'),.05))

    def test_each_field_endpoint_uses_its_own_geometry_backend(self):
        seen=[]
        def sample(cp,points):
            seen.append(os.environ['LAGFORCE_SEEDED_POWER_CELLS']);return {'density':np.ones(len(points))}
        with patch.dict(os.environ,{}),patch('compare_physical_budget.sampled_fields',side_effect=sample):
            values=sampled_pair([{},{}],[dict(power_cell_evaluation_backend='active'),dict(power_cell_evaluation_backend='seeded_active')],np.array([[.2,.3]]))
        self.assertEqual(seen,['0','1']);self.assertEqual(len(values),2)

    def test_physical_scales_and_distinct_density_goal(self):
        scales,goals=scales_and_goals(dict(problem='steady_vortex',amplitude=.5,physical_goal=.01,physical_goals=dict(density_field=.001)))
        self.assertEqual(scales['density'],1.);self.assertEqual(scales['velocity'],.5)
        self.assertEqual(goals,dict(density=.001,velocity=.01,pressure=.01))
        scales,_=scales_and_goals(dict(amplitude=.02,physical_goal=.01));self.assertEqual(scales,dict(density=.02,velocity=.02,pressure=.02))

    def test_model_and_exact_initial_array_guard(self):
        cfg=dict(n=4,amplitude=.02,kinetic_smoothing_coefficient=1.)
        data=dict(t=0.,mass=[1.],state=[2.],factors=[[[0.,0.],[0.,0.]]])
        compatible([cfg,cfg],[data,data])
        with self.assertRaises(ValueError):compatible([cfg,dict(cfg,kinetic_smoothing_coefficient=0.)],[data,data])
        with self.assertRaises(ValueError):compatible([cfg,cfg],[data,dict(data,mass=[1.00001])])

    def test_vector_and_scalar_l1_normalization(self):
        left=dict(density=np.array([1.,3.]),velocity=np.array([[3.,4.],[0.,0.]]),pressure=np.array([1.,1.]))
        right={k:np.zeros_like(v) for k,v in left.items()}
        self.assertEqual(field_difference(left,right,dict(density=2.,velocity=.5,pressure=1.)),dict(density=1.,velocity=5.,pressure=1.))


if __name__=='__main__':unittest.main()

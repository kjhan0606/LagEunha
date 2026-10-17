import unittest
import numpy as np
from common_contact_time_refinement import digest,difference,uniform_grid,observed_orders,end_time_agrees
from common_contact_fixed_time_comparison import validate_cohort_time
from common_contact_long_time import install_input


class TimeComparisonTests(unittest.TestCase):
    def test_endpoint_roundoff_is_distinct_from_a_different_final_time(self):
        self.assertTrue(end_time_agrees(.0005099999999999997,.00051,.0004))
        self.assertFalse(end_time_agrees(.00051+1e-12,.00051,.0004))
        self.assertFalse(end_time_agrees(float('nan'),.00051,.0004))

    def test_long_comparison_resumes_only_identical_input(self):
        import tempfile
        from pathlib import Path
        original=dict(config=dict(n=3),t=.00011,state=[1.],mass=[1.],factors=[[[0.]]])
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'inputs/coarse/checkpoint.json'
            self.assertTrue(install_input(path,original))
            self.assertFalse(install_input(path,original))
            with self.assertRaises(RuntimeError):install_input(path,dict(original,state=[2.]))

    def test_explicit_cohort_time_remains_enforced(self):
        validate_cohort_time(.00011,.00011)
        with self.assertRaises(RuntimeError):validate_cohort_time(1e-5,.00011)
        with self.assertRaises(RuntimeError):validate_cohort_time(.00011,float('nan'))

    def test_rejection_free_run_may_still_be_nonuniform(self):
        cp=dict(steps=4,rejected=0)
        self.assertTrue(uniform_grid(cp,[.25]*4,1.,4))
        self.assertFalse(uniform_grid(cp,[.25,.2,.25,.3],1.,4))
        self.assertFalse(uniform_grid(dict(steps=5,rejected=0),[.2]*5,1.,4))
        self.assertFalse(uniform_grid(dict(steps=4,rejected=1),[.25]*4,1.,4))

    def test_order_requires_uniform_grid_and_resolved_differences(self):
        d=[dict(x=4e-6,w=4e-15),dict(x=1e-6,w=1e-15)]
        self.assertEqual(observed_orders(d,1e-13,True),dict(x=2.,w=None))
        self.assertEqual(observed_orders(d,1e-13,False),dict(x=None,w=None))

    def test_gauge_weight_and_mass_normalization(self):
        mass=np.array([1.,2.]);a=np.zeros((2,20));b=a.copy()
        b[:,2]=3.;b[:,3]=mass*1e-4
        ca=dict(state=np.r_[a.ravel(),np.zeros(4)].tolist())
        cb=dict(state=np.r_[b.ravel(),np.zeros(4)].tolist())
        d=difference(ca,cb,mass)
        self.assertEqual(d['gauge_fixed_weight'],0.)
        self.assertAlmostEqual(d['raw_moments_over_mass'],1e-4)
        self.assertEqual(d['generator'],0.)

    def test_state_digest_survives_json_roundtrip(self):
        import json
        state=np.array([.1,1/3,-.25])
        self.assertEqual(digest(state),digest(json.loads(json.dumps(state.tolist()))))
        self.assertNotEqual(digest(state),digest(state+1e-10))


if __name__=='__main__':unittest.main()

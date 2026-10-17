import unittest
import numpy as np
from audit_3d_history import summarize,validate_chain


class HistoryAuditTests(unittest.TestCase):
    def test_accumulate_cellwise_and_use_full_duration(self):
        rows=[dict(signed=np.array([1.,0.]),absolute=np.array([1.,0.]),quad_difference=np.zeros(2),relative_rate=2.,generator=.01),
              dict(signed=np.array([-1.,2.]),absolute=np.array([1.,2.]),quad_difference=np.ones(2)*.1,relative_rate=1.,generator=.02)]
        d=summarize(rows,np.ones(2),10.)
        self.assertEqual(d['max_absolute_relative_mass'],2.)
        self.assertEqual(d['duration_scaled_peak_rate'],20.)
        self.assertEqual(d['maximum_generator_path_defect'],.02)
        self.assertEqual(d['summed_nested_quadrature_difference'],.1)

    def test_missing_time_and_changed_mass_rejected(self):
        state=dict(sites=np.zeros((1,3)),weights=np.zeros(1),raw=np.ones((1,35)),factors=np.eye(3)[None])
        record=dict(state,t=0.);end={k:v.copy() for k,v in state.items()}
        validate_chain(record,end,state,0.,np.ones(1))
        with self.assertRaisesRegex(ValueError,'gap'):validate_chain(record,end,state,1.,np.ones(1))
        end['raw'][0,0]=2.
        with self.assertRaises(AssertionError):validate_chain(record,end,state,0.,np.ones(1))

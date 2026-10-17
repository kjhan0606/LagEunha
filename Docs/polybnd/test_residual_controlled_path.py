import unittest
import numpy as np
from residual_controlled_path import integrate
from replay_dense_weight_path import segment_value


class ControllerTests(unittest.TestCase):
    def test_rejection_propagates_accepted_endpoint(self):
        segments,end,attempts=integrate(lambda t,y:np.ones_like(y),[2.],0.,1.,
            lambda s:s.h,.08)
        self.assertTrue(any(not r['accepted'] for r in attempts))
        self.assertTrue(all(s.h<=.08 for s in segments))
        np.testing.assert_allclose(end,[3.],rtol=0,atol=2e-14)
        for left,right in zip(segments,segments[1:]):
            self.assertAlmostEqual(left.t_old+left.h,right.t_old)
            np.testing.assert_array_equal(segment_value(left,right.t_old),right.y_old)

    def test_exhausted_refinement_fails(self):
        with self.assertRaisesRegex(RuntimeError,'depth exhausted'):
            integrate(lambda t,y:np.zeros_like(y),[1.],0.,1.,lambda s:1.,.1,max_depth=1)

    def test_nonfinite_defect_fails(self):
        with self.assertRaisesRegex(RuntimeError,'nonfinite'):
            integrate(lambda t,y:np.zeros_like(y),[1.],0.,1.,lambda s:np.nan,.1)


if __name__=='__main__':unittest.main()

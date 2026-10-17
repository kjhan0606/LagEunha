import unittest
import numpy as np
from scipy.integrate import solve_ivp
from continue_residual_path import unpack, continuity
from local_weight_defect import serialize
from replay_dense_weight_path import segment_value


class ContinuationTests(unittest.TestCase):
    def test_serialized_continuation(self):
        fun=lambda t,y:np.array([1.,-2.])
        left=solve_ivp(fun,(0.,.2),[.7,1.],method='DOP853',dense_output=True,rtol=1e-12,atol=1e-14)
        restored=unpack(dict(segments=[serialize(s) for s in left.sol.interpolants]))
        end=segment_value(restored[-1],.2)
        right=solve_ivp(fun,(.2,.4),end,method='DOP853',dense_output=True,rtol=1e-12,atol=1e-14)
        report=continuity(restored+right.sol.interpolants,2)
        self.assertLess(report['max_state_jump'],1e-14)
        self.assertLess(report['max_weight_derivative_jump'],1e-12)
        np.testing.assert_allclose(right.y[:,-1],[1.1,.2],atol=1e-14,rtol=0)

    def test_gap_rejected(self):
        sol=solve_ivp(lambda t,y:y*0,(0.,1.),[1.],method='DOP853',dense_output=True)
        segments=unpack(dict(segments=[serialize(s) for s in sol.sol.interpolants]))
        segments[1].t_old+=.01
        with self.assertRaises(AssertionError):continuity(segments,1)


if __name__=='__main__':unittest.main()

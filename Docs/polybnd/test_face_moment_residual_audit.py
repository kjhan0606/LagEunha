import unittest
import numpy as np
from face_moment_residual_audit import minimax


class ResidualTests(unittest.TestCase):
    def test_known_relaxation_and_dual_bound(self):
        p=dict(k=1,A=np.array([[1.]]),b=np.array([.5]),bounds=[(.1,1.01),(0,None)])
        result=minimax(p,np.array([[1.]]),np.array([.8]))
        self.assertTrue(result['gate'])
        self.assertAlmostEqual(result['minimum_normalized_residual'],.3,places=12)
        self.assertAlmostEqual(result['dual_objective'],.3,places=12)


if __name__=='__main__':unittest.main()

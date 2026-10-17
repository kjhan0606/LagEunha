import unittest
import numpy as np
from bounded_moment_feasibility import minimum_moment_residual,box_dual_lower_bound


class BoundedMomentFeasibilityTests(unittest.TestCase):
    def test_repaired_dual_never_exceeds_known_optimum(self):
        matrix=np.array([[1.,0.,-1.],[0.,1.,-1.],[-1.,0.,-1.],[0.,-1.,-1.]])
        rhs=np.array([1.2,-.1,-1.2,.1]);cost=np.array([0.,0.,1.])
        d=box_dual_lower_bound(matrix,rhs,cost,[(0.,1.),(0.,1.),(0.,None)],np.array([1e6,-2e6,3e6,-4e6]))
        self.assertGreaterEqual(d['slack_reduced_cost'],0.)
        self.assertLessEqual(d['lower_bound'],.2)

    def test_feasible_box_has_validated_witness(self):
        r=minimum_moment_residual(np.eye(2),np.array([.3,.7]),[(0.,1.)]*2)
        self.assertEqual(r['classification'],'moment_witness')
        self.assertLess(r['actual_moment_residual'],1e-10)

    def test_incompatible_box_has_positive_dual_certificate(self):
        r=minimum_moment_residual(np.eye(2),np.array([1.2,-.1]),[(0.,1.)]*2)
        self.assertEqual(r['classification'],'finite_basis_infeasible_certificate')
        self.assertAlmostEqual(r['dual_lower_bound'],.2)
        self.assertLess(r['duality_gap'],1e-10)


if __name__=='__main__':unittest.main()

import unittest
import numpy as np
from density_failure_diagnosis import diagnose_density_failure
from moment_aligned_density_audit import polygon_moments


class DensityFailureDiagnosisTests(unittest.TestCase):
    def geometry(self):return polygon_moments(np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]]))

    def problem(self,b,bounds=None):
        return dict(low=0.,high=1.,raw=np.asarray([self.geometry()]),k=2,A=np.eye(2),
            b=np.asarray(b),bounds=[(0.,1.)]*2 if bounds is None else bounds)

    def test_bound_conflict_is_independent_of_basis(self):
        g=self.geometry();half=polygon_moments(np.array([[.5,0.],[1.,0.],[1.,1.],[.5,1.]]))
        p=self.problem([.3,.7]);p.update(low=1.,high=50.5)
        d=diagnose_density_failure(p,g+99*half)
        self.assertEqual(d['classification'],'density_bounds_incompatible_with_moments')
        self.assertGreater(d['density_bound_certificate']['required_upper_at_least'],90.)

    def test_finite_basis_has_positive_residual_certificate(self):
        d=diagnose_density_failure(self.problem([1.2,-.1]),.5*self.geometry())
        self.assertEqual(d['classification'],'finite_basis_infeasible_certificate')
        self.assertGreater(d['minimum_residual_diagnosis']['dual_lower_bound'],.19)

    def test_gate_witness_does_not_claim_exact_solver_feasibility(self):
        d=diagnose_density_failure(self.problem([.3,.7]),.5*self.geometry())
        self.assertEqual(d['classification'],'moment_gate_witness_solver_unresolved')

    def test_fixed_face_bounds_are_part_of_the_finite_problem(self):
        d=diagnose_density_failure(self.problem([.3,.7],[(0.,0.),(1.,1.)]),.5*self.geometry(),'face')
        self.assertEqual(d['classification'],'finite_basis_and_trace_infeasible_certificate')


if __name__=='__main__':unittest.main()

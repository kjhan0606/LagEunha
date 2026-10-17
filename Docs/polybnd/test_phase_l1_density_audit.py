import unittest
import numpy as np
from subcell_density_audit import triangle_moments
from phase_l1_density_audit import phase_fit


class PhaseObjectiveTests(unittest.TestCase):
    def test_objective_matches_integral_with_fixed_density(self):
        v=np.array([[0.,0.],[1.,0.],[0.,1.]])
        callback=lambda parts:np.full(len(parts),1.5)
        callback.phase_data=lambda parts:(1.,3.,np.full(len(parts),.25))
        result,_,rho=phase_fit(v,2*triangle_moments(v),1,None,(2.,2.),callback)
        self.assertTrue(result['feasible'] and result['moment_gate'] and result['bound_gate'])
        # Constant rho=2: integral |2-rho_prior|=area; divide by mass=2*area.
        self.assertAlmostEqual(result['inferred_prior_L1_relative'],.5,places=12)
        self.assertLess(result['objective_error'],1e-12)
        np.testing.assert_allclose(rho,2.,atol=1e-12)

    def test_constant_prior_zero_objective(self):
        v=np.array([[0.,0.],[1.,0.],[0.,1.]])
        callback=lambda parts:np.full(len(parts),2.)
        callback.phase_data=lambda parts:(2.,2.,np.zeros(len(parts)))
        result,_,_=phase_fit(v,2*triangle_moments(v),1,None,(1.,3.),callback)
        self.assertTrue(result['feasible'] and result['moment_gate'])
        self.assertAlmostEqual(result['inferred_prior_L1_relative'],0.,places=12)


if __name__=='__main__':unittest.main()

import unittest
import numpy as np
from subcell_density_audit import fit,triangle_moments
from bounded_density_audit import minimum_relaxation


class BoundTests(unittest.TestCase):
    def setUp(self):
        self.v=np.array([[0.,0.],[1.,0.],[0.,1.]])
        self.target=2*triangle_moments(self.v)

    def test_constant_at_bounds(self):
        result,_,rho=fit(self.v,self.target,1,lambda x:np.full(len(x),2.),(2.,2.))
        self.assertTrue(result['feasible'] and result['moment_gate'] and result['bound_gate'])
        np.testing.assert_allclose(rho,2.,atol=1e-12)

    def test_upper_bound_conflicts_with_mass(self):
        result,_,rho=fit(self.v,self.target,1,lambda x:np.ones(len(x)),(.5,1.))
        self.assertFalse(result['feasible']);self.assertIsNone(rho)

    def test_minimum_relaxation_from_mass(self):
        result=minimum_relaxation(self.v,self.target,1,1.,1.)
        self.assertTrue(result['feasible'] and result['gate'])
        self.assertAlmostEqual(result['minimum_relaxation'],1.,places=10)


if __name__=='__main__':unittest.main()

import unittest
import numpy as np
from constrained_weight_path import PathSolution

class ConstrainedPathTests(unittest.TestCase):
    def test_stage_mass_and_rate(self):
        p=PathSolution(.3,732);a=p.at(.5)
        self.assertLess(a['mass_error'],1e-11)
        self.assertLess(a['instantaneous_mass_flux'],1e-13)
        self.assertGreater(a['max_face'],1e-5)
    def test_constraint_derivative(self):
        p=PathSolution(.3,731);t=.5;h=1e-5
        fd=(p.at(t+h)['w']-p.at(t-h)['w'])/(2*h)
        np.testing.assert_allclose(fd,p.at(t)['z'],atol=2e-9,rtol=0)
if __name__=='__main__':unittest.main()

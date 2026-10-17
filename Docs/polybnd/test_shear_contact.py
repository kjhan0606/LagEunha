import unittest
from adaptive_moment_batch import setup
from adaptive_coupled_mass import configuration
from shear_contact_check import affine_errors


class ShearContactTests(unittest.TestCase):
    def test_initial_affine_identities(self):
        for ratio in (10.,100.):
            mass,y=setup('shear_contact',4,ratio)
            self.assertLess(max(affine_errors(y,mass).values()),1e-12)

    def test_detect_momentum_error(self):
        mass,y=setup('shear_contact',3,100.)
        y[9]+=mass[0]*.01
        self.assertGreater(affine_errors(y,mass)['p'],.009)

    def test_invalid_ratio(self):
        with self.assertRaises(ValueError):configuration('shear_contact',3,0.,.01,5e-11)

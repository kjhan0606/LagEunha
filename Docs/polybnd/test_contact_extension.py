import unittest
import numpy as np
from adaptive_coupled_mass import configuration
from adaptive_moment_batch import setup
from kinetic_factor_rhs import rhs


class ContactTests(unittest.TestCase):
    def test_configuration(self):
        self.assertEqual(configuration('contact',4,0.,.05,5e-11)['kind'],'contact')

    def test_stationary_contact_rhs(self):
        mass,y=setup('contact',4,0.)
        factors=[np.zeros((2,2)) for _ in mass]
        dy,d=rhs(0.,y,mass,covariance_factors=factors)
        self.assertLess(float(np.max(abs(dy))),1e-10)
        self.assertLess(d['mass_flux'],1e-12)
        self.assertAlmostEqual(d['min_pressure'],.6,places=10)

import unittest
from audit_kepler_energy_step import balance


class EnergyBalanceTests(unittest.TestCase):
    def test_gravity_work_cancels_potential_change(self):
        self.assertEqual(balance(2.,3.,-4.,-5.,0.,0.),0.)

    def test_incoming_boundary_exchange_is_subtracted(self):
        self.assertEqual(balance(2.,5.,-4.,-3.,2.,2.),0.)
        self.assertEqual(balance(2.,5.,-4.,-3.,2.,1.),1.)

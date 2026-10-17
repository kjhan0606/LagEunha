import unittest
import numpy as np
from extension_time_convergence import compare


class ConvergenceTests(unittest.TestCase):
    def test_second_order(self):
        reference = np.zeros((9, 20))
        result = compare([np.ones_like(reference)/k**2 for k in (16, 32, 64)], reference)
        for row in result.values():
            self.assertAlmostEqual(row['observed_order'], 2.)
            self.assertTrue(row['second_order_gate'])

    def test_unresolved_is_not_a_pass(self):
        state = np.zeros((9, 20))
        for row in compare([state, state, state], state).values():
            self.assertFalse(row['resolved'])
            self.assertIsNone(row['second_order_gate'])

    def test_first_order_rejected(self):
        reference = np.zeros((9, 20))
        for row in compare([np.ones_like(reference)/k for k in (16, 32, 64)], reference).values():
            self.assertFalse(row['second_order_gate'])

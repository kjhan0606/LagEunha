import unittest
import numpy as np
from audit_saved_dense_3d import integrate


class SavedDenseAuditTests(unittest.TestCase):
    def test_cubic_integral_and_nested_rule(self):
        s=np.linspace(0,1,17)
        values=np.column_stack([1+s+s*s+s**3, -2*s**3])
        expected=.2*np.array([1+.5+1/3+.25,-.5])
        np.testing.assert_allclose(integrate(values,.2),expected,atol=1e-15,rtol=0)
        np.testing.assert_allclose(integrate(values[::2],.2),expected,atol=1e-15,rtol=0)

    def test_invalid_counts(self):
        for n in (0,1,2,4):
            with self.assertRaises(ValueError):integrate(np.zeros((n,2)),.2)

import unittest
from types import SimpleNamespace
import numpy as np
from scipy.integrate._ivp.rk import Dop853DenseOutput
from replay_dense_weight_path import segment_value


class ReplayTests(unittest.TestCase):
    def test_values_match_scipy(self):
        rng = np.random.default_rng(418289)
        for a, b in ((0., 1.), (.3, .30001)):
            old, coefficients = rng.normal(size=13), rng.normal(size=(7, 13))
            reference = Dop853DenseOutput(a, b, old, coefficients)
            saved = SimpleNamespace(t_old=a, h=b-a, y_old=old, F=coefficients)
            for t in np.linspace(a, b, 19):
                np.testing.assert_array_equal(segment_value(saved, t), reference(t))


if __name__ == '__main__':
    unittest.main()

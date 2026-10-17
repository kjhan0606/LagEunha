import unittest
from finalize_shear_time_isolation import sensitivity


class SensitivityTests(unittest.TestCase):
    def test_relative_change(self):
        rows=[dict(interface=dict(cell_average_density_L1_relative=v)) for v in (2e-5,1.1e-5,1.01e-5)]
        result=sensitivity(rows,dict(cell_average_density_L1_relative=1e-5))
        self.assertAlmostEqual(result['finest_change_relative_to_bias'],.01)

    def test_zero_baseline_rejected(self):
        with self.assertRaises(ValueError):sensitivity([],dict(cell_average_density_L1_relative=0.))

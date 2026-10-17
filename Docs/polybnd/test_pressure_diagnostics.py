import unittest
import numpy as np
from simple_wave_evolution import pressure_diagnostics,update_pressure_peak


class PressureDiagnosticsTests(unittest.TestCase):
    def test_only_small_scalar_reconstruction_metadata_is_saved(self):
        d=pressure_diagnostics(dict(pressure_rates=np.ones(10),pressure_polynomial_limited_cells=2,pressure_reconstruction='polynomial_cell_mean_degree4',min_pressure=.5))
        self.assertNotIn('pressure_rates',d);self.assertNotIn('min_pressure',d);self.assertEqual(d['pressure_polynomial_limited_cells'],2)

    def test_peaks_are_accumulated_without_changing_existing_diagnostics(self):
        peak=dict(generator=1e-8);update_pressure_peak(peak,[dict(pressure_polynomial_limited_cells=2,pressure_polynomial_min_positivity_scale=.8,pressure_polynomial_neighbor_count_deficient_cells=1)])
        update_pressure_peak(peak,[dict(pressure_polynomial_limited_cells=1,pressure_polynomial_min_positivity_scale=.9)])
        self.assertEqual(peak['pressure_limited_cells_max'],2);self.assertEqual(peak['pressure_positivity_scale_min'],.8);self.assertEqual(peak['pressure_diagnostics_steps_observed'],2);self.assertEqual(peak['generator'],1e-8)


if __name__=='__main__':unittest.main()

import unittest
from unittest.mock import patch
import json
import numpy as np
from active_power_cells import power_cells_active
from weighted_mass_patch import power_cells


class ActivePowerCellTests(unittest.TestCase):
    def test_seeded_clipping_keeps_all_nonlocal_weighted_constraints(self):
        rng=np.random.default_rng(20261017);domain=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        for count in [4,16,64]:
            sites=rng.uniform(.02,.98,(count,2));weights=rng.normal(0.,.12,count)
            for seed_count in [2,12]:
                original=power_cells_active(sites,weights,domain);seeded=power_cells_active(sites,weights,domain,seed_count)
                for a,b in zip(original,seeded):self.assertLess(a.symmetric_difference(b).area,3e-12)

    def test_seeded_geometry_requires_explicit_backend_configuration(self):
        from simple_wave_evolution import run_with_executor
        env=dict(LAGFORCE_VECTOR_KINETIC_ROOT='0',LAGFORCE_QUADRATIC_VELOCITY='0',LAGFORCE_VECTOR_QUADRATIC_BUDGET='0',LAGFORCE_QUADRATIC_PRESSURE='0',LAGFORCE_VECTOR_SMOOTH_PATH='0',LAGFORCE_VECTOR_MOMENT_FLUX='0',LAGFORCE_ACTIVE_POWER_CELLS='1',LAGFORCE_SEEDED_POWER_CELLS='1')
        with patch.dict('os.environ',env),patch('simple_wave_evolution.Path.read_bytes',return_value=json.dumps({}).encode()):
            with self.assertRaisesRegex(ValueError,'power cell evaluation backend'):run_with_executor()

    def test_seeded_regular_grid_preserves_cells_area_and_weight_gauge(self):
        x,y=np.meshgrid((np.arange(7)+.5)/7,(np.arange(7)+.5)/7);sites=np.column_stack([x.ravel(),y.ravel()]);weights=np.zeros(len(sites));domain=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        a=power_cells_active(sites,weights,domain,12);b=power_cells_active(sites,weights+12.,domain,12)
        for p,q in zip(a,b):self.assertAlmostEqual(p.area,1/49,places=13);self.assertLess(p.symmetric_difference(q).area,2e-13)

    def test_all_constraints_match_exhaustive_weighted_clipping(self):
        rng=np.random.default_rng(20261016);domain=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        for count in (4,16,64):
            sites=rng.uniform(.02,.98,(count,2))
            for weights in (np.zeros(count),rng.uniform(-.002,.002,count)):
                original=power_cells(sites,weights,domain);active=power_cells_active(sites,weights,domain)
                self.assertLess(abs(sum(p.area for p in active)-1.),1e-12)
                for a,b in zip(original,active):self.assertLess(a.symmetric_difference(b).area,2e-12)

    def test_weight_gauge_and_empty_cell(self):
        domain=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        sites=np.array([[.25,.25],[.75,.25],[.25,.75],[.75,.75]])
        weights=np.array([-2.,0.,0.,0.])
        a=power_cells_active(sites,weights,domain);b=power_cells_active(sites,weights+100.,domain)
        self.assertTrue(a[0].is_empty)
        for p,q in zip(a,b):self.assertLess(p.symmetric_difference(q).area,1e-12)

    def test_duplicate_or_nonfinite_sites_are_rejected(self):
        domain=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        with self.assertRaises(ValueError):power_cells_active(np.zeros((2,2)),np.zeros(2),domain)
        with self.assertRaises(ValueError):power_cells_active(np.array([[np.nan,0.]]),np.zeros(1),domain)


if __name__=='__main__':unittest.main()

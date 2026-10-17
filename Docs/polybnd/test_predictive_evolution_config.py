import unittest
import numpy as np
from adaptive_moment_batch import setup
from common_interface_contact import initialize
import predictive_evolution_batch as batch
from predictive_boundary_density import BoundaryDensity
from predictive_mass_path import evaluate_path
import kinetic_factor_rhs as kinetic


class EvolutionConfigurationTests(unittest.TestCase):
    def test_explicit_grid_selection_keeps_legacy_default(self):
        import json,tempfile
        from pathlib import Path
        from unittest.mock import patch
        from common_contact_evolution import selected_cases,CASES
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'experiment_config.json'
            with patch('common_contact_evolution.Path',return_value=p):
                self.assertEqual(selected_cases(),CASES)
                c=dict(name='n8_angle23',n=8,angle=23.)
                p.write_text(json.dumps(dict(cases=[c])))
                self.assertEqual(selected_cases(),[('n8_angle23',8,23.)])
                p.write_text(json.dumps(dict(cases=[c,c])))
                with self.assertRaises(ValueError):selected_cases()

    def config(self):
        return dict(kind='common_planar_contact',n=3,amplitude=100.,angle=23.,smooth_phase_floor=True,
            density_provider=dict(mode='moment',segments=16,radial=1e-7,core_level=2,align_core=True,
                construct_prior_witness=True,lp_retry=True))

    def test_common_mass_dispatch_does_not_use_cellwise_mass(self):
        cfg=self.config();mass,y=batch.initialize_case(cfg);expected,state,_,_=initialize(3,100.,23.)
        np.testing.assert_array_equal(mass,expected);np.testing.assert_array_equal(y,state)
        old,_=setup('shear_contact',3,100.)
        self.assertGreater(float(np.max(abs(old-mass))),1.)

    def test_legacy_dispatch_and_provider_defaults(self):
        cfg=dict(kind='shear_contact',n=3,amplitude=10.)
        mass,y=batch.initialize_case(cfg);expected,state=setup('shear_contact',3,10.)
        np.testing.assert_array_equal(mass,expected);np.testing.assert_array_equal(y,state)
        self.assertEqual(batch.density_options(cfg,False),dict(mode='neighbor',segments=64,radial=1e-7,
            vector_geometry=False,smooth_phase_floor=False))
        with self.assertRaises(ValueError):batch.initialize_case(dict(cfg,kind='unrecognized'))

    def test_two_phase_reference_on_many_cells(self):
        import moment_ale_evolve as model
        from common_contact_evolution import reference_moments
        mass,y,refs,phases=initialize(3,100.,23.);a=y[:-4].reshape(-1,20)
        cells,_=model.geometry(a[:,:2],a[:,2],0.)
        exact,coverage,_=reference_moments(cells,refs,phases,0.)
        np.testing.assert_allclose(exact,np.column_stack([mass,a[:,3:9]]),rtol=1e-12,atol=1e-13)
        self.assertLess(coverage,1e-14)

    def test_path_worker_uses_same_density_as_stage(self):
        cfg=self.config();mass,y=batch.initialize_case(cfg);options=batch.density_options(cfg,True)
        provider=BoundaryDensity(**options)
        batch.initialize_path_worker(3,100.,True,True,options)
        worker=batch._path_provider
        for key,value in options.items():self.assertEqual(getattr(worker,key),value)
        L=[np.zeros((2,2)) for _ in mass]
        dy,diag=kinetic.rhs(0.,y,mass,covariance_factors=L,face_density_provider=provider,
            face_segments=provider.segments,face_order=4,exact_bin_transport=True)
        rates=dy[:-4].reshape(-1,20)
        residual,_=evaluate_path(0.,y,L,mass,(rates[:,:2],rates[:,2]),worker)
        np.testing.assert_allclose(residual,diag['mass_rates'],rtol=0,atol=1e-11)
        self.assertLess(float(np.max(abs(residual)/mass)),1e-10)


if __name__=='__main__':unittest.main()

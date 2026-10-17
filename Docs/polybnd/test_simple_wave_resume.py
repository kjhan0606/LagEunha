import hashlib,json,tempfile,unittest
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from simple_wave_evolution import finish_pending_commit,preserved_initial


class SimpleWaveResumeTests(unittest.TestCase):
    def test_saved_approval_completes_commit_without_rhs(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory);state=np.array([.1,.2,.3]);identity='example'
            cp=dict(config_sha256=identity,t=.001,state=[.4,.5,.6])
            trial=dict(record=dict(t=0.,y=state),commit=dict(checkpoint=cp,accepted=True,
                archive='step_1_trial.json',sample_archive='step_1_samples.json'))
            atomic(p/'pending.json',trial);atomic(p/'pending_samples.json',dict(sample=1))
            atomic(p/'source_rhs.json',dict(state_sha256=hashlib.sha256(state.tobytes()).hexdigest()))
            # Recover both before and after the atomic checkpoint write.
            atomic(p/'checkpoint.json',cp)
            finish_pending_commit(p/'pending.json',p/'checkpoint.json',p/'pending_samples.json',identity)
            self.assertEqual(json.loads((p/'checkpoint.json').read_text()),cp)
            self.assertTrue((p/'step_1_trial.json').exists());self.assertTrue((p/'step_1_samples.json').exists())
            self.assertFalse((p/'source_rhs.json').exists());self.assertFalse((p/'pending.json').exists())
            finish_pending_commit(p/'pending.json',p/'checkpoint.json',p/'pending_samples.json',identity)

    def test_wrong_identity_cannot_commit(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory);atomic(p/'pending.json',dict(commit=dict(checkpoint=dict(config_sha256='wrong'))))
            with self.assertRaises(ValueError):finish_pending_commit(p/'pending.json',p/'checkpoint.json',p/'samples.json','expected')
            self.assertFalse((p/'checkpoint.json').exists())

    def test_preserved_physical_source_identity_and_mass(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'initial.json';source=dict(t=0.,mass=[1.],state=[0.]*24,factors=[[[0.,0.],[0.,0.]]],
                physical_initial_condition=dict(kind='nonlinear_simple_wave',n=1,amplitude=.02))
            atomic(p,source);cfg=dict(n=1,amplitude=.02,end=.05,preserved_initial_state=str(p),
                preserved_initial_state_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
            saved,mass,state,L=preserved_initial(cfg);np.testing.assert_array_equal(mass,[1.])
            with self.assertRaises(ValueError):preserved_initial(dict(cfg,amplitude=.01))

    def test_stationary_vortex_source_cannot_be_mislabeled_as_wave(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'initial.json';source=dict(t=0.,mass=[1.],state=[0.]*24,factors=[[[0.,0.],[0.,0.]]],
                physical_initial_condition=dict(kind='steady_vortex',n=1,amplitude=.5))
            atomic(p,source);cfg=dict(n=1,amplitude=.5,end=.05,problem='steady_vortex',preserved_initial_state=str(p),
                preserved_initial_state_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
            saved,mass,state,L=preserved_initial(cfg);np.testing.assert_array_equal(mass,[1.])
            with self.assertRaises(ValueError):preserved_initial(dict(cfg,problem='nonlinear_simple_wave'))

    def test_saved_fluid_stages_complete_only_the_geometry(self):
        from simple_wave_evolution import complete_geometry_trial
        from unittest.mock import patch
        y=np.arange(24,dtype=float);end=y+1.;fake=y.copy();fake[:3]=[.2,.3,.4]
        trial=dict(record=dict(t=0.,dt=.001,y=y.copy(),L=[np.eye(2)],d0=np.zeros(24),dm=np.zeros(24),diag={}),
            new=end.copy(),newL=[np.eye(2)],stage_mass=1e-18,geometry_pending=True)
        with patch('geometry_mass_lift.integrate',return_value=dict(segments=[])),patch('geometry_mass_lift.lifted_state',return_value=(fake,[np.eye(2)],np.zeros(24))):
            result=complete_geometry_trial(trial,np.array([1.]),{})
        np.testing.assert_array_equal(result['new'][3:],end[3:]);np.testing.assert_array_equal(result['new'][:3],fake[:3])
        np.testing.assert_array_equal(result['record']['y'],y);self.assertNotIn('geometry_pending',result)
        self.assertEqual(result['stage_mass'],1e-18)


if __name__=='__main__':unittest.main()

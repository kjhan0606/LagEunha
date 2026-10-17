import hashlib,json,tempfile,unittest
from pathlib import Path
from copy import deepcopy
import fcntl
import numpy as np
from adaptive_coupled_mass import atomic
from resume_numerical_acceleration import prepare,sha,reidentify


class NumericalAccelerationResumeTests(unittest.TestCase):
    def test_additional_equivalent_backends_preserve_primary_state(self):
        with tempfile.TemporaryDirectory() as folder:
            source,target,control=self.setup_source(folder);before=json.loads((source/'checkpoint.json').read_text());backends=dict(mass_path_evaluation_backend='cell_batched',moment_flux_evaluation_backend='batched',power_cell_evaluation_backend='seeded_active');control.update(backends);atomic(target/'numerical_resume_config.json',control)
            marker=prepare(target);cfg=json.loads((target/'simple_wave_config.json').read_text());cp=json.loads((target/'checkpoint.json').read_text())
            for k,v in backends.items():self.assertEqual(cfg[k],v);self.assertEqual(marker['numerical_backend_changes'][k],v)
            for k,v in before.items():
                if k!='config_sha256':self.assertEqual(cp[k],v)

    def test_unknown_acceleration_backend_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            _,target,control=self.setup_source(folder);control['mass_path_evaluation_backend']='unverified';atomic(target/'numerical_resume_config.json',control)
            with self.assertRaisesRegex(ValueError,'unsupported acceleration backend'):prepare(target)

    def test_physical_covariance_equivalence_is_checked_when_provided(self):
        with tempfile.TemporaryDirectory() as folder:
            _,target,control=self.setup_source(folder);proof=Path(control['equivalence_report']);data=json.loads(proof.read_text());data['physical_covariance_rate_maximum_absolute_difference']=1e-3;atomic(proof,data);control['equivalence_report_sha256']=sha(proof);atomic(target/'numerical_resume_config.json',control)
            with self.assertRaisesRegex(ValueError,'physical covariance'):prepare(target)

    def setup_source(self,folder):
        p=Path(folder);source=p/'source';target=p/'target';source.mkdir();target.mkdir();(source/'OWNER').write_text('LagEunha fake source fixture\n');(source/'execution.txt').write_text('exit=143\n');(source/'.run.lock').touch()
        cfg=dict(n=1,amplitude=.02,end=.2,geometry_path='quadratic_midpoint',velocity_reconstruction='quadratic_budget_projected',path_mass_budget=1e-5);atomic(source/'simple_wave_config.json',cfg);identity=sha(source/'simple_wave_config.json')
        state=np.arange(24,dtype=float).tolist();initial=dict(t=0.,mass=[1.],state=state,factors=[np.eye(2).tolist()],config_sha256=identity);cp=dict(initial,t=.1,next_dt=.01,steps=2,rejected=1,absolute=[1e-7],signed=[-1e-8],quad_delta=[1e-11],peak=dict(generator=1e-8))
        atomic(source/'initial_state.json',initial);atomic(source/'checkpoint.json',cp)
        atomic(source/'pending.json',dict(record=dict(t=.1,dt=.01,y=state,d0=state,dm=state),new=state,newL=initial['factors']))
        atomic(source/'pending_samples.json',{'0':dict(flux=[0.],diag=dict(generator_defect=0.))})
        proof=p/'proof.json';atomic(proof,dict(completed=True,rhs_maximum_absolute_difference=1e-18,independent_mass_flux_maximum_absolute_difference=1e-19))
        control=dict(source=str(source),source_config_sha256=identity,source_initial_sha256=sha(source/'initial_state.json'),source_checkpoint_sha256=sha(source/'checkpoint.json'),equivalence_report=str(proof),equivalence_report_sha256=sha(proof));atomic(target/'numerical_resume_config.json',control)
        return source,target,control

    def test_primary_arrays_accumulators_and_both_pending_files_are_preserved(self):
        with tempfile.TemporaryDirectory() as folder:
            source,target,_=self.setup_source(folder);before={n:sha(source/n) for n in ['initial_state.json','checkpoint.json','pending.json','pending_samples.json']};old=json.loads((source/'checkpoint.json').read_text());r=prepare(target)
            current=json.loads((target/'checkpoint.json').read_text())
            for k in old:
                if k!='config_sha256':self.assertEqual(current[k],old[k])
            self.assertEqual(current['config_sha256'],sha(target/'simple_wave_config.json'));self.assertTrue(r['pending_fluid_stages_reused']);self.assertTrue(r['pending_sample_cache_reused'])
            for n,h in before.items():self.assertEqual(sha(source/n),h)
            self.assertEqual(sha(target/'pending.json'),before['pending.json']);self.assertEqual(sha(target/'pending_samples.json'),before['pending_samples.json'])

    def test_idempotent_preparation_never_restores_an_older_approved_state(self):
        with tempfile.TemporaryDirectory() as folder:
            _,target,_=self.setup_source(folder);prepare(target);p=target/'checkpoint.json';cp=json.loads(p.read_text());cp['state'][0]+=.001;cp['t']+=.001;cp['steps']+=1;atomic(p,cp);h=sha(p)
            self.assertTrue(prepare(target)['already_prepared']);self.assertEqual(sha(p),h)

    def test_source_hash_guard_and_running_lock(self):
        with tempfile.TemporaryDirectory() as folder:
            source,target,control=self.setup_source(folder)
            with (source/'.run.lock').open('a') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                with self.assertRaisesRegex(ValueError,'still locked'):prepare(target)
            control['source_checkpoint_sha256']='0'*64;atomic(target/'numerical_resume_config.json',control)
            with self.assertRaisesRegex(ValueError,'hash mismatch'):prepare(target)

    def test_saved_approved_commit_is_finished_before_identity_change(self):
        with tempfile.TemporaryDirectory() as folder:
            source,target,control=self.setup_source(folder);cp=json.loads((source/'checkpoint.json').read_text());new=deepcopy(cp);new['state'][0]+=.01;new.update(t=.11,steps=3)
            trial=json.loads((source/'pending.json').read_text());trial['new']=new['state'];trial['commit']=dict(checkpoint=new,accepted=True,archive='step_3_trial.json',sample_archive='step_3_samples.json');atomic(source/'pending.json',trial)
            r=prepare(target);current=json.loads((target/'checkpoint.json').read_text());self.assertEqual(current['state'],new['state']);self.assertEqual(current['steps'],3);self.assertTrue(r['source_pending_commit_completed']);self.assertFalse((target/'pending.json').exists());self.assertTrue((target/'step_3_trial.json').exists())

    def test_reidentify_only_changes_the_configuration_key(self):
        original=dict(config_sha256='old',state=[1.,2.],mass=[1.]);changed=reidentify(original,'new');self.assertEqual(original['config_sha256'],'old');self.assertEqual(changed['state'],original['state']);self.assertEqual(changed['config_sha256'],'new')

    def test_cached_rhs_from_another_state_is_refused(self):
        with tempfile.TemporaryDirectory() as folder:
            source,target,_=self.setup_source(folder);atomic(source/'source_rhs.json',dict(t=.1,state_sha256='0'*64))
            with self.assertRaisesRegex(ValueError,'cached source RHS'):prepare(target)

    def test_horner_backend_updates_only_a_higher_pressure_configuration(self):
        with tempfile.TemporaryDirectory() as folder:
            source,target,control=self.setup_source(folder);cfg=json.loads((source/'simple_wave_config.json').read_text());cfg['pressure_reconstruction']='polynomial_cell_mean_degree4';atomic(source/'simple_wave_config.json',cfg);identity=sha(source/'simple_wave_config.json')
            for name in ['initial_state.json','checkpoint.json']:
                d=json.loads((source/name).read_text());atomic(source/name,reidentify(d,identity))
            control.update(source_config_sha256=identity,source_initial_sha256=sha(source/'initial_state.json'),source_checkpoint_sha256=sha(source/'checkpoint.json'),pressure_evaluation_backend='horner');atomic(target/'numerical_resume_config.json',control)
            prepare(target);self.assertEqual(json.loads((target/'simple_wave_config.json').read_text())['pressure_evaluation_backend'],'horner')


if __name__=='__main__':unittest.main()

"""Policies and dispatch failure paths; no actual Slurm command or GPU work."""
from contextlib import contextmanager
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import time
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import run_kepler_checkpoint_validation_20261029 as runner


class CheckpointDispatchTests(unittest.TestCase):
    def environment(self):
        return {'SLURM_JOB_ID':'438957','SLURM_CPUS_PER_TASK':'2',
                'SLURM_MEM_PER_NODE':'32768','SLURM_JOB_NUM_NODES':'1',
                'SLURM_JOB_END_TIME':str(int(time.time())+10000)}

    def gates(self):
        return ({'time':3.,'checks':{str(i):True for i in range(15)}},
                {'time':3.,'checks':{str(i):True for i in range(9)}})

    def test_gate_requires_all_checks_exact_time_and_resources(self):
        g,r=self.gates();env=self.environment()
        runner.gate_policy(g,r,env,'438957')
        for bad_g,bad_r,bad_env in [({**g,'checks':{}},r,env),
                                    ({**g,'time':2.99},r,env),
                                    (g,{**r,'checks':{'a':True}},env),
                                    (g,r,{**env,'SLURM_MEM_PER_NODE':'8192'}),
                                    (g,r,{**env,'SLURM_JOB_ID':'439115'})]:
            with self.assertRaises(ValueError):runner.gate_policy(bad_g,bad_r,bad_env,'438957')

    def test_step_cleanup_scope_excludes_parent_batch_and_other_jobs(self):
        self.assertEqual(runner.numeric_steps('438957.batch|R\n438957.extern|R\n438957.7|R','438957'),['438957.7'])
        for line in ['438957|R','999999.7|R','438957.some_unknown_step|R']:
            with self.assertRaises(ValueError):runner.numeric_steps(line,'438957')
        with patch.object(runner,'steps',side_effect=[['438957.7'],[]]),patch.object(runner,'command') as cmd,patch.object(runner.time,'sleep') as sleep:
            clean,_=runner.cleanup_steps('438957')
            self.assertTrue(clean)
            cmd.assert_called_once_with('scancel','--signal=TERM','438957.7')
            sleep.assert_called_once_with(60)

    def test_uses_step_ID_format_supported_by_current_Slurm(self):
        with patch.object(runner,'command',return_value='438957.7\n438957.batch\n438957.extern\n') as cmd:
            self.assertEqual(runner.steps('438957'),['438957.7'])
            cmd.assert_called_once_with('squeue','--steps','--jobs','438957','--noheader','--format=%i')

    @contextmanager
    def fixture(self):
        previous=Path.cwd()
        with TemporaryDirectory() as directory:
            base=Path(directory);stage=base/'stage';root=base/'benchmark';stage.mkdir();root.mkdir()
            g,r=self.gates()
            (stage/'audit_gate.json').write_text(json.dumps(g));(stage/'continuation_gate.json').write_text(json.dumps(r))
            (stage/'stage.json').write_text(json.dumps({'next_job':'/unchanged/next.slurm'}))
            config={'stage':str(stage),'benchmark':str(root),'parent_job':438957,'queued_job':439115,
                    'expected_next_job':'/unchanged/next.slurm','cases':['control40'],'frozen_hashes':{}}
            path=base/'config.json';path.write_text(json.dumps(config))
            os.chdir(stage)
            try:
                with patch.dict(os.environ,self.environment()):yield stage,root,path
            finally:os.chdir(previous)

    def test_already_running_standalone_never_held_or_cancelled(self):
        with self.fixture() as (stage,root,path):
            commands=[]
            def command(*args):
                commands.append(args)
                if args[:4]==('scontrol','show','job','438957'):return 'JobState=RUNNING'
                if args[:4]==('scontrol','show','job','439115'):return 'JobState=RUNNING RunTime=00:00:10'
                raise AssertionError(args)
            with patch.object(runner,'command',side_effect=command),patch.object(runner,'steps',return_value=[]),patch.object(runner.subprocess,'Popen') as launch:
                self.assertEqual(runner.run(path),0);launch.assert_not_called()
            self.assertFalse(any(c[0]=='scancel' for c in commands))
            self.assertFalse(any(c[:2]==('scontrol','hold') for c in commands))
            self.assertEqual(json.loads((stage/'stage.json').read_text())['next_job'],'/unchanged/next.slurm')

    def test_unconfirmed_cancellation_cannot_launch_inline(self):
        with self.fixture() as (stage,root,path):
            commands=[]
            def command(*args):
                commands.append(args)
                if args[:4]==('scontrol','show','job','438957'):return 'JobState=RUNNING'
                if args[:4]==('scontrol','show','job','439115'):
                    return f'JobState=PENDING RunTime=00:00:00 UserId=test({os.getuid()}) Reason=QOSMaxGRESPerUser WorkDir={root} JobName=KepVertexHaloPair'
                return ''
            with patch.object(runner,'command',side_effect=command),patch.object(runner,'steps',return_value=[]),patch.object(runner.subprocess,'Popen') as launch:
                self.assertEqual(runner.run(path),0);launch.assert_not_called()
            self.assertIn(('scancel','--ctld','--state=PENDING','439115'),commands)
            self.assertIn(('scontrol','release','439115'),commands)
            report=json.loads((stage/'post_checkpoint_validation.json').read_text())
            self.assertFalse(report['queued_request_cancelled']);self.assertTrue(report['continuation_metadata_unchanged'])

    def test_duplicate_dispatch_is_not_retried(self):
        with self.fixture() as (stage,root,path):
            (stage/'post_checkpoint_validation.lock').write_text('already claimed')
            with patch.object(runner,'command') as cmd:
                self.assertEqual(runner.run(path),0);cmd.assert_not_called()

    def test_deadline_uses_actual_allocation_end(self):
        self.assertEqual(runner.remaining_seconds({'SLURM_JOB_END_TIME':'5000'},{},4000),1000)
        with self.assertRaises(ValueError):runner.remaining_seconds({}, {},4000)

    def child_result(self, code, clean):
        with self.fixture() as (stage,root,path):
            queued_queries=0
            def command(*args):
                nonlocal queued_queries
                if args[:4]==('scontrol','show','job','438957'):return 'JobState=RUNNING'
                if args[:4]==('scontrol','show','job','439115'):
                    queued_queries+=1
                    state='CANCELLED' if queued_queries>=3 else 'PENDING'
                    return f'JobState={state} RunTime=00:00:00 UserId=test({os.getuid()}) Reason=QOSMaxGRESPerUser WorkDir={root} JobName=KepVertexHaloPair'
                return ''
            process=Mock(pid=999999)
            process.wait.return_value=code
            process.poll.return_value=code
            def launch(*args,**kwargs):
                if code==0:(root/'comparison.json').write_text(json.dumps({'passed':True}))
                return process
            preflight=SimpleNamespace(stdout=b'checks passed',stderr=b'',check_returncode=lambda:None)
            with (patch.object(runner,'command',side_effect=command),
                  patch.object(runner,'steps',return_value=[]),
                  patch.object(runner.subprocess,'run',return_value=preflight),
                  patch.object(runner.subprocess,'Popen',side_effect=launch),
                  patch.object(runner,'cleanup_steps',return_value=(clean,[]))):
                result=runner.run(path)
            report=json.loads((stage/'post_checkpoint_validation.json').read_text())
            self.assertTrue(report['continuation_metadata_unchanged'])
            self.assertFalse(report['automatic_promotion'])
            self.assertTrue(report['queued_request_cancelled'])
            return result,report

    def test_completed_auxiliary_preserves_original_continuation(self):
        result,report=self.child_result(0,True)
        self.assertEqual(result,0);self.assertTrue(report['auxiliary_passed'])
        result,report=self.child_result(92,True)
        self.assertEqual(result,0);self.assertFalse(report['auxiliary_passed'])

    def test_unreaped_auxiliary_steps_cannot_advance_pipeline(self):
        result,report=self.child_result(92,False)
        self.assertEqual(result,32);self.assertFalse(report['auxiliary_passed'])


if __name__=='__main__':unittest.main()

"""Mocked orchestration contracts: never launches GPU, MPI, Slurm or real signals."""
from contextlib import ExitStack
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import run_kepler_halo40_promotion_20261029 as engine
import run_kepler_checkpoint_validation_20261029 as guard
import audit_kepler_padding_reuse_pair_20261028 as common


class PromotionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.stage, self.root, self.chain, self.cache = [self.base/n for n in ['stage','root','chain','cache']]
        for p in [self.stage,self.root,self.chain,self.cache]:p.mkdir()
        self.original = str(self.base/'old/job.slurm')
        self.candidate = str(self.chain/'to4p0/job.slurm')
        (self.stage/'stage.json').write_text(json.dumps({'next_job':self.original,'target_time':3.}))
        self.before = (self.stage/'stage.json').read_bytes()
        self.checkpoint = self.base/'khout.023456.dat'
        self.checkpoint.write_bytes(b'accepted-checkpoint')
        (self.stage/'audit_gate.json').write_text(json.dumps({'time':3.,'endpoint':str(self.checkpoint),
            'checks':{str(i):True for i in range(15)}}))
        (self.stage/'continuation_gate.json').write_text(json.dumps({'time':3.,'checks':{str(i):True for i in range(9)}}))
        (self.stage/'post_checkpoint_validation.json').write_text(json.dumps({'status':'AUXILIARY_FAILED_OR_INCOMPLETE',
            'auxiliary_process_group':123,'auxiliary_steps_reaped':True}))
        self.gpu = {'host':'syn103','GPU_UUID':'GPU-mock','GPU_model':'A100 mock'}
        self.env = {'SLURM_JOB_ID':'438957','SLURM_NODELIST':'syn103','CUDA_VISIBLE_DEVICES':'0',
                    'SLURM_CPUS_PER_TASK':'2','SLURM_MEM_PER_NODE':'32768','SLURM_JOB_NUM_NODES':'1'}
        for name in ['canonical.exe','candidate.exe','params.dat','initial.dat',
                     'audit_kepler_polar_long_gate_20261027.py','audit_kepler_native_radial_20261026.py']:
            (self.root/name).write_bytes(b'mock-input')
        (self.root/'templates').mkdir()
        template = '\n'.join(['export '+k+'='+v for k,v in {
            'GFS_KEPLER_HALO_CELLS':'40','GFS_GPU_VERTEX_HALO_CERT':'0','GFS_REUSE_TESS_PADDING':'0',
            'GFS_GPU_LOCAL_PERIODIC_SYNC':'0','GFS_GEOM_FORCE_DUMP_MAXCALL':'4','HYDRO_TSTOP':'2.66',
            'EUNHA_DUMP_DT':'.05','GFS_EXACT_STOP':'1'}.items()])+'''
"$PY" /mock/record_kepler_vertex_case_context_20261029.py
mpirun ./eunha2.exe params.dat 23008 > log
"$PY" ./audit_kepler_polar_long_gate_20261027.py ./khout.000000.dat "$PWD" --time 2.66 --output audit_endpoint_gate.json
'''
        for name in ['canonical.sh','control.sh']:(self.root/'templates'/name).write_text(template)
        self.cases = ['fresh_reference40','combined40','old_a','combined_b','combined_b2','old_a2']
        for name in self.cases:(self.root/name).mkdir()
        (self.cache/'launch_context.json').write_text(json.dumps({**self.env,
            **dict(zip(engine.FLAGS,map(str,(0,40,0,0))))}))
        (self.cache/'audit_endpoint_gate.json').write_text(json.dumps({'checks':{'mock':True}}))
        (self.cache/'eunha2.exe').write_bytes(b'candidate')
        for i in range(1,5):(self.cache/f'geom_force_rank000_call{i:03d}.bin').write_bytes(b'mock')
        self.config = {'parent_job':'438957','stage':str(self.stage),'root':str(self.root),
            'candidate_root':str(self.chain),'expected_next_job':self.original,
            'cached_control':str(self.cache),'candidate_binary_sha256':engine.sha(self.cache/'eunha2.exe'),
            'candidate_template':{},'static_cases':self.cases,'frozen_hashes':{}}
        self.config_path=self.base/'config.json'
        self.write_config()

    def write_config(self):self.config_path.write_text(json.dumps(self.config))

    def fake_load(self,name):
        bridge=name.startswith('bridge_')
        metrics={'steps':1 if bridge else 59,'time':3.001 if bridge else 2.66,
                 'wall_seconds':75. if name.startswith('combined') else 100.}
        return None,{'steps':[('3.001','23457','.001')]},metrics

    def fake_process(self,*args,**kwargs):
        case=Path(kwargs['cwd']);name=case.name;test=self
        class Process:
            pid=99999999
            def poll(self):return 1
            def wait(self,timeout=None):
                if test.mode=='timeout':raise subprocess.TimeoutExpired('mock',timeout)
                if test.mode=='native_failure':return 1
                endpoint=case/'endpoint.dat';endpoint.write_bytes(b'endpoint')
                (case/'audit_endpoint_gate.json').write_text(json.dumps({'endpoint':str(endpoint)}))
                combo=name.startswith('combined') or name=='bridge_combo'
                flags=(0,40,1,1) if combo else (0,40,0,0)
                context={**test.gpu,'SLURM_JOB_ID':'438957','OMP_NUM_THREADS':'2',
                         'GFS_EXACT_STOP':'0' if name.startswith('bridge') else '1','GFS_DIAG_DT_MAX':'',
                         **dict(zip(engine.FLAGS,map(str,flags)))}
                (case/'launch_context.json').write_text(json.dumps(context))
                (case/'log').write_text('[TESS_PADDING_REUSE] rank=0 x reuse=1 refresh_validation=1\n'
                    '[LOCAL_PERIODIC_SYNC] rank=0 device alias refresh\n' if combo else '')
                return 0
        return Process()

    def invoke(self,mode='success',cleanup=True,stale=False,residual=False):
        self.mode=mode
        old_cwd=Path.cwd();os.chdir(self.stage);self.addCleanup(os.chdir,old_cwd)
        actual_atomic=engine.atomic
        def write(path,value):
            if stale and path.name=='net_wall_comparison.json':value={**value,'passed':False}
            return actual_atomic(path,value)
        with ExitStack() as stack:
            stack.enter_context(patch.dict(os.environ,self.env,clear=True))
            stack.enter_context(patch.object(engine,'device',return_value=self.gpu))
            stack.enter_context(patch.object(engine,'snapshot',side_effect=lambda p:
                {'time':3. if Path(p)==self.checkpoint else 3.001,'dt':.001,'dtold':.001}))
            stack.enter_context(patch.object(engine,'atomic',side_effect=write))
            stack.enter_context(patch.object(guard,'steps',return_value=['438957.99'] if residual else []))
            stack.enter_context(patch.object(guard,'command',return_value=f'JobState=RUNNING UserId=mock({os.getuid()})'))
            stack.enter_context(patch.object(guard,'remaining_seconds',return_value=7200))
            self.cleaner=stack.enter_context(patch.object(guard,'cleanup_steps',return_value=(cleanup,[])))
            stack.enter_context(patch.object(guard,'stop_group',return_value=-15))
            stack.enter_context(patch.object(subprocess,'check_output',return_value=''))
            stack.enter_context(patch.object(subprocess,'run',return_value=subprocess.CompletedProcess([],0)))
            self.launcher=stack.enter_context(patch.object(subprocess,'Popen',side_effect=self.fake_process))
            stack.enter_context(patch.object(common,'compare',return_value={'passed':True}))
            stack.enter_context(patch.object(common,'load',side_effect=self.fake_load))
            stack.enter_context(patch('audit_kepler_geometry_force_dump_20261029.compare_dumps',return_value={'passed':True}))
            stack.enter_context(patch('audit_kepler_restart_bytes_20261029.compare',return_value={'passed':True}))
            stack.enter_context(patch('kepler_halo40_candidate_chain_20261029.bind_chain',return_value='b'*64))
            return engine.run(self.config_path)

    def test_success_commits_only_after_all_gates(self):
        self.assertEqual(self.invoke(),0)
        self.assertEqual(json.loads((self.stage/'stage.json').read_text())['next_job'],self.candidate)
        self.assertTrue(json.loads((self.stage/'halo40_post_checkpoint.json').read_text())['promoted'])
        self.assertEqual(self.launcher.call_count,8)  # Reused reference, 1 combo, 4 ABBA, 3 bridge.

    def test_timeout_retains_original_after_cleanup(self):
        self.assertEqual(self.invoke('timeout'),0)
        self.assertEqual((self.stage/'stage.json').read_bytes(),self.before)
        self.cleaner.assert_called_once()

    def test_unreaped_steps_stop_at_checkpoint(self):
        self.assertEqual(self.invoke('native_failure',cleanup=False),32)
        self.assertEqual((self.stage/'stage.json').read_bytes(),self.before)

    def test_original_residual_rank_is_not_signalled_or_advanced(self):
        self.assertEqual(self.invoke(residual=True),32)
        self.launcher.assert_not_called();self.cleaner.assert_not_called()

    def test_partial_output_is_not_replayed(self):
        (self.root/'combined40/log').write_text('partial')
        self.assertEqual(self.invoke(),0)
        self.launcher.assert_not_called()
        self.assertEqual((self.stage/'stage.json').read_bytes(),self.before)

    def test_stale_failed_receipt_prevents_promotion(self):
        self.assertEqual(self.invoke(stale=True),0)
        self.assertEqual((self.stage/'stage.json').read_bytes(),self.before)

    def test_duplicate_live_and_stale_reports_stop(self):
        h=engine.sha(self.config_path)
        good={'finished_at':'now','cleanup_verified':True,'safe_to_advance':True,
              'config_sha256':h,'status':'ORIGINAL_NEXT_RETAINED'}
        self.assertEqual(engine.duplicate_result({},h,self.original,self.original,self.candidate),32)
        self.assertEqual(engine.duplicate_result(good,h,self.original,self.original,self.candidate),0)
        self.assertEqual(engine.duplicate_result({**good,'config_sha256':'wrong'},h,self.original,self.original,self.candidate),32)
        self.assertEqual(engine.duplicate_result({**good,'status':'OPERATIONAL_STOP'},h,self.original,self.original,self.candidate),32)

    def test_compare_and_swap_preserves_operator_change(self):
        (self.stage/'stage.json').write_text(json.dumps({'next_job':'operator_choice'}))
        eligibility={key:'a'*64 for key in ['halo40','performance','bridge','checkpoint_sha256']}
        with self.assertRaises(ValueError):engine.promotion_transaction(self.stage,self.before,Path(self.candidate),eligibility,'b'*64)
        self.assertEqual(json.loads((self.stage/'stage.json').read_text())['next_job'],'operator_choice')

    def test_bridge_keeps_natural_step_no_cap(self):
        engine.make_case(self.root,'natural','off',self.checkpoint,23456,3.0000000000000004,4,bridge=True)
        text=(self.root/'natural/run.sh').read_text()
        self.assertIn('GFS_EXACT_STOP=0',text)
        self.assertIn('unset GFS_DIAG_DT_MAX',text)
        self.assertIn('GFS_DIAG_CHECKPOINT_STEP=23457',text)
        self.assertIn('"$PWD/audit_view" --time "$actual_time"',text)


if __name__=='__main__':unittest.main()

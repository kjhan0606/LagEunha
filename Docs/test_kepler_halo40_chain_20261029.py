"""Candidate graph/input binding tests; no native execution or scheduler access."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import kepler_halo40_candidate_chain_20261029 as chain


class CandidateChainTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.base=Path(self.tmp.name);self.old=self.base/'old';self.old.mkdir()
        self.new=self.base/'candidate';self.binary=self.base/'binary'
        self.binary.write_bytes(b'candidate')
        for helper in ['prepare_restart.py','audit_kepler_polar_long_gate_20261027.py',
                       'audit_kepler_native_radial_20261026.py',
                       'check_radial_and_next.before_checkpoint_validation_20261029.py']:
            (self.old/helper).write_text('# original gate fixture\n')
        previous=3.;previous_name='to3p0'
        for i,(name,target) in enumerate(chain.STAGES):
            p=self.old/name;p.mkdir()
            next_name=chain.STAGES[i+1][0] if i+1<len(chain.STAGES) else None
            meta={'initial_time':previous,'target_time':target,
                  'parent_gate':str(self.old/previous_name/'audit_gate.json'),
                  'next_job':str(self.old/next_name/'job.slurm') if next_name else None,
                  'checkpoint':None,'checkpoint_step':None,'checkpoint_sha256':None}
            (p/'stage.json').write_text(json.dumps(meta))
            (p/'params.dat').write_text('original physics parameters\n')
            (p/'khout.000000.dat').write_bytes(b'initial')
            (p/'job.slurm').write_text('''#!/usr/bin/env bash
#SBATCH --nodes=1
PY=/home/kjhan/miniconda3/bin/python3
export GFS_GPU_LOCAL_PERIODIC_SYNC=0
TIMEFORMAT='wall=%R'
''')
            previous,previous_name=target,name

    def prepare(self):
        with patch.object(chain,'CANDIDATE_SHA',chain.sha(self.binary)),patch.object(chain,'BASE_CHECKER_SHA',chain.sha(
            self.old/'check_radial_and_next.before_checkpoint_validation_20261029.py')):
            return chain.stage_chain(self.old,self.new,self.binary)

    def test_original_preserved_and_candidate_stops_inputs_exact(self):
        originals={name:(self.old/name/'stage.json').read_bytes() for name,_ in chain.STAGES}
        manifest=self.prepare();chain.verify_files(self.new,manifest)
        for name,target in chain.STAGES:
            self.assertEqual((self.old/name/'stage.json').read_bytes(),originals[name])
            meta=json.loads((self.new/name/'stage.json').read_text())
            self.assertEqual(meta['target_time'],target)
            self.assertEqual((self.new/name/'params.dat').read_bytes(),(self.old/name/'params.dat').read_bytes())
            text=(self.new/name/'job.slurm').read_text()
            self.assertIn('GFS_GPU_LOCAL_PERIODIC_SYNC=1',text)
            self.assertIn('GFS_REUSE_TESS_PADDING=1 GFS_GPU_VERTEX_HALO_CERT=0',text)
            self.assertIn('--nodelist=syn103',text)
        first=json.loads((self.new/'to4p0/stage.json').read_text())
        self.assertEqual(first['parent_gate'],str(self.old/'to3p0/audit_gate.json'))
        self.assertFalse((self.new/'chain_files.sha256').exists())

    def test_bind_uses_original_checkpoint_then_rejects_mutation(self):
        manifest=self.prepare();checkpoint=self.base/'accepted3.dat';checkpoint.write_bytes(b'original accepted3')
        result=chain.bind_chain(self.new,manifest,checkpoint,12345,chain.sha(checkpoint),
            {'host':'syn103','GPU_UUID':'GPU-mock','GPU_model':'A100 mock'})
        self.assertEqual(len(result),64)
        first=json.loads((self.new/'to4p0/stage.json').read_text())
        self.assertEqual(first['checkpoint'],str(checkpoint))
        self.assertEqual(first['checkpoint_step'],12345)
        self.assertTrue(json.loads((self.new/'device_binding.json').read_text())['bound'])
        bound=json.loads((self.new/'bound_manifest.json').read_text())['files']
        chain.verify_files(self.new,bound)
        (self.new/'to5p0/params.dat').write_text('changed')
        with self.assertRaises(ValueError):chain.verify_files(self.new,bound)

    def test_existing_candidate_output_prevents_binding(self):
        manifest=self.prepare();checkpoint=self.base/'accepted3.dat';checkpoint.write_bytes(b'checkpoint')
        (self.new/'to4p0/log').write_text('existing run')
        with self.assertRaises(ValueError):chain.bind_chain(self.new,manifest,checkpoint,12345,chain.sha(checkpoint),{})


if __name__=='__main__':unittest.main()

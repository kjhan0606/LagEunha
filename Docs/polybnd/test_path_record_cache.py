import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock,patch
import numpy as np
import predictive_evolution_batch as batch
from predictive_path_refinement import assess_adaptive


class PathRecordCacheTests(unittest.TestCase):
    def record(self,dt=.1):
        return dict(t=0.,dt=dt,y=[1.],L=[[[0.,0.],[0.,0.]]],d0=[2.],dm=[3.],
            diag=dict(cov_A=[[[0.,0.],[0.,0.]]],cov_Qfactor=[[[0.],[0.]]]))

    def test_same_record_is_loaded_once_per_worker(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'pending.json';record=self.record()
            path.write_text(json.dumps(dict(record=record)));expected=batch.record_digest(record)
            with patch.multiple(batch,_path_record_key=None,_path_record_data=None),patch.object(batch,'path_sample',return_value=7) as sample:
                self.assertEqual(batch.path_sample_file(path,expected,np.ones(1),0),7)
                path.unlink()
                self.assertEqual(batch.path_sample_file(path,expected,np.ones(1),32),7)
                self.assertEqual(sample.call_count,2)
                np.testing.assert_array_equal(sample.call_args.args[0]['y'],record['y'])

    def test_a_new_digest_loads_the_new_trial(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'pending.json'
            with patch.multiple(batch,_path_record_key=None,_path_record_data=None),patch.object(batch,'path_sample',return_value=0) as sample:
                for dt in (.1,.2):
                    record=self.record(dt);path.write_text(json.dumps(dict(record=record)))
                    batch.path_sample_file(path,batch.record_digest(record),np.ones(1),0)
                    self.assertEqual(sample.call_args.args[0]['dt'],dt)

    def test_wrong_digest_is_rejected_before_evaluation(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'pending.json';path.write_text(json.dumps(dict(record=self.record())))
            with patch.multiple(batch,_path_record_key=None,_path_record_data=None),patch.object(batch,'path_sample') as sample:
                with self.assertRaises(RuntimeError):batch.path_sample_file(path,'wrong',np.ones(1),0)
                sample.assert_not_called()

    def test_submission_sends_path_and_digest(self):
        executor=Mock();executor.submit.side_effect=[object(),object()];record=self.record()
        batch.submit_path_samples(executor,record,np.ones(1),[2,4],2.,'pending.json')
        self.assertEqual(executor.submit.call_count,2)
        for call,index in zip(executor.submit.call_args_list,(1.,2.)):
            args=call.args
            self.assertIs(args[0],batch.path_sample_file)
            self.assertEqual(args[2],batch.record_digest(record));self.assertEqual(args[4],index)
            self.assertEqual(len(args),5);self.assertIsInstance(args[1],str);self.assertIsInstance(args[2],str)

    def test_larger_lattice_reuses_old_noncoarse_samples(self):
        def value(s):return dict(mass_rate=[1e-8*s**3],generator=0.,pressure=.6,factor_R=0.)
        record=dict(t=0.,dt=1e-6);original={str(i):value(i/64) for i in range(65)};seen=[]
        def sample(rec,mass,index,provider=None):seen.append(index);return value(index/64)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'cache.json';old=value(1/16384)
            path.write_text(json.dumps(dict(version=2,record_sha256=batch.record_digest(record),samples={'1':old})))
            with patch.object(batch,'path_sample',side_effect=sample):
                result=assess_adaptive(record,np.ones(1),1e-4,None,original,path,'unit',min_nodes=2049,max_nodes=2049)
            saved=json.loads(path.read_text())
            self.assertEqual(result['refined_nodes'],2049)
            self.assertEqual(saved['lattice'],32768);self.assertEqual(saved['samples']['2'],old)
            self.assertNotIn(1/256,seen)
            self.assertAlmostEqual(result['signed'][0],2.5e-15,places=27)

    def test_unsupported_lattice_is_rejected(self):
        with self.assertRaises(ValueError):assess_adaptive({},np.ones(1),1.,None,{},'unused','unit',max_nodes=512)


if __name__=='__main__':unittest.main()

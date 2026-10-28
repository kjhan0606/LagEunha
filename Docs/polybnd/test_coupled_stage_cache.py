import pickle
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from coupled_stage_cache import StageCache
import potential_moment_coupled_2d as coupled


class CoupledStageCacheTests(unittest.TestCase):
    def assert_exact(self, expected, actual):
        self.assertIs(type(expected),type(actual))
        if isinstance(expected,np.ndarray):
            self.assertEqual(expected.dtype,actual.dtype)
            self.assertEqual(expected.shape,actual.shape)
            self.assertEqual(expected.tobytes(),actual.tobytes())
        elif isinstance(expected,dict):
            self.assertEqual(list(expected),list(actual))
            for key in expected:self.assert_exact(expected[key],actual[key])
        elif isinstance(expected,(list,tuple)):
            self.assertEqual(len(expected),len(actual))
            for a,b in zip(expected,actual):self.assert_exact(a,b)
        else:self.assertEqual(pickle.dumps(expected,protocol=5),pickle.dumps(actual,protocol=5))

    def test_exact_types_and_corruption_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            cache=StageCache(tmp,{'source':'a','dt':.01})
            value=(np.arange(6,dtype=np.longdouble).reshape(2,3),[np.eye(2)],{'flag':True})
            cache('gravity_first',lambda:value)
            restored=cache('gravity_first',lambda:self.fail('recomputed'))
            self.assertEqual(pickle.dumps(value,protocol=5),pickle.dumps(restored,protocol=5))
            with self.assertRaisesRegex(ValueError,'identity'):
                StageCache(tmp,{'source':'b','dt':.01})('gravity_first',lambda:None)
            p=Path(tmp)/'gravity_first.cache';data=p.read_bytes();p.write_bytes(data[:-1]+bytes([data[-1]^1]))
            with self.assertRaisesRegex(ValueError,'checksum'):cache('gravity_first',lambda:None)

    def test_failed_stage_not_committed_and_partial_file_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            cache=StageCache(tmp,{'source':'a'})
            def interrupted():raise RuntimeError('interrupted')
            with self.assertRaisesRegex(RuntimeError,'interrupted'):cache('hydro_first',interrupted)
            self.assertFalse((Path(tmp)/'hydro_first.cache').exists())
            (Path(tmp)/'hydro_first.tmp').write_bytes(b'partial')
            self.assertEqual(cache('hydro_first',lambda:42),42)

    def test_stage_name_is_bound(self):
        with tempfile.TemporaryDirectory() as tmp:
            cache=StageCache(tmp,{'source':'a'});cache('gravity_first',lambda:1)
            (Path(tmp)/'gravity_last.cache').write_bytes((Path(tmp)/'gravity_first.cache').read_bytes())
            with self.assertRaisesRegex(ValueError,'identity'):cache('gravity_last',lambda:2)

    def test_interrupted_coupled_step_reuses_only_completed_stages(self):
        # Controlled numerical kernels isolate persistence from the hydro model.
        y=np.arange(24,dtype=float)/100;mass=np.ones(1);L=np.eye(2)[None];U=np.array([-.5])
        calls=[]
        def kick(t,y,m,L,g,h,*args):
            calls.append(('kick',t));return y+h,L.copy(),{'gravity_work':h}
        def hydro(t,y,m,L,U,*args):
            calls.append(('hydro',t))
            diag=dict(generator_velocity_defect=0.,mass_flux=0.,potential_boundary_rate=.1,gravity_work_rate=.2,
                cov_linearization_factors=L,cov_A=L,cov_K=L,cov_Qfactor=L)
            return y*.01,U*.02,diag
        def compute(cache=None):
            return coupled.step(0.,y,mass,L,U,.01,None,None,None,return_record=True,stage_cache=cache)
        with patch.object(coupled,'apply',side_effect=kick),patch.object(coupled,'stage',side_effect=hydro),patch.object(coupled,'mask_at',return_value=np.array([False])),patch.object(coupled,'restore_frozen'),patch.object(coupled,'rebuild',side_effect=lambda y,m,L:y),patch.object(coupled,'advance_factors_sylvester',side_effect=lambda L,d,h:L.copy()):
            reference=compute()
            for stop in ('hydro_first','hydro_middle','gravity_last'):
                with self.subTest(stop=stop),tempfile.TemporaryDirectory() as tmp:
                    cache=StageCache(tmp,{'source':'fixture'})
                    def interrupt(name,fn):
                        if name==stop:raise RuntimeError('interrupted')
                        return cache(name,fn)
                    with self.assertRaisesRegex(RuntimeError,'interrupted'):compute(interrupt)
                    calls.clear();actual=compute(StageCache(tmp,{'source':'fixture'}))
                    self.assert_exact(reference,actual)
                    self.assertNotIn(('kick',0.),calls)
                    if stop!='hydro_first':self.assertNotIn(('hydro',0.),calls)
                    if stop=='gravity_last':self.assertNotIn(('hydro',.005),calls)


if __name__=='__main__':unittest.main()

import os,hashlib,unittest
from unittest.mock import patch
import numpy as np
from resume_saved_branch import resume_task


class SavedBranchResumeTests(unittest.TestCase):
    def fixtures(self):
        source=dict(t=1e-5,state=[1.,2.],mass=[3.,4.])
        cfg=dict(case='trial',interval='short',start=1e-5,end=1.1e-4,duration=1e-4,initial_steps=4,
            source_state_sha256=hashlib.sha256(np.asarray(source['state']).tobytes()).hexdigest())
        return source,dict(config=cfg,t=2e-5,mass=[3.,4.])

    def test_task_preserves_original_interval_after_checkpoint_resume(self):
        source,cp=self.fixtures()
        with patch.dict(os.environ,{'LAGFORCE_ANALYTIC_CAP':'0','LAGFORCE_CENTERED_MOMENT_FIT':'0','LAGFORCE_ACTIVE_POWER_CELLS':'0'}):
            self.assertEqual(resume_task('trial',source,cp),('trial','short',1e-4,4))
            cp['mass'][0]+=1
            with self.assertRaises(ValueError):resume_task('trial',source,cp)

    def test_mixed_numerical_flags_or_changed_source_are_rejected(self):
        source,cp=self.fixtures()
        with patch.dict(os.environ,{'LAGFORCE_ANALYTIC_CAP':'1'}):
            with self.assertRaises(ValueError):resume_task('trial',source,cp)
        source['state'][0]+=1
        with self.assertRaises(ValueError):resume_task('trial',source,cp)


if __name__=='__main__':unittest.main()

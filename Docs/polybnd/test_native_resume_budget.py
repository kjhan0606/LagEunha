import unittest
import numpy as np
from resume_native_physical_budget import native_trial


class NativeResumeBudgetTests(unittest.TestCase):
    def test_reuses_fluid_stages_and_changes_only_geometry(self):
        y=np.arange(24,dtype=float);dm=np.ones(24);new=y+2
        trial=dict(record=dict(y=y.copy(),dm=dm.copy(),dt=.001,geometry_lift=dict(example='old')),new=new.copy(),newL=[np.eye(2)],stage_mass=1e-18,geometry_pending=True)
        result=native_trial(trial)
        np.testing.assert_array_equal(result['new'][3:],new[3:]);np.testing.assert_array_equal(result['new'][:3],y[:3]+.001)
        np.testing.assert_array_equal(result['record']['y'],y);np.testing.assert_array_equal(result['record']['dm'],dm)
        self.assertNotIn('geometry_lift',result['record']);self.assertNotIn('geometry_pending',result);self.assertEqual(result['stage_mass'],1e-18)


if __name__=='__main__':unittest.main()

import copy
import unittest
import numpy as np
from compare_temporal_states import compare


class TemporalStatesTests(unittest.TestCase):
    def fixtures(self):
        return [dict(t=.25,mass=[1.],state=(np.arange(24.)+h*h).tolist(),
                     potential_moments=[-1.+h*h],factors=[np.eye(2).tolist()])
                for h in (.1,.05,.025)]

    def test_second_order_and_zero_differences(self):
        result=compare(self.fixtures(),[.1,.05,.025])
        for key,row in result['variables'].items():
            if key=='covariance':
                self.assertFalse(row['differences_resolved'])
                self.assertIsNone(row['observed_order'])
            else:
                self.assertAlmostEqual(row['observed_order'],2.,places=9)
                self.assertAlmostEqual(row['difference_alignment'],1.)
        self.assertFalse(result['full_validation'])

    def test_mismatches_rejected(self):
        for key,value in [('t',.3),('mass',[2.]),('state',[0.]),('potential_moments',[float('nan')])]:
            states=copy.deepcopy(self.fixtures());states[1][key]=value
            with self.subTest(key=key),self.assertRaises((ValueError,AssertionError)):
                compare(states,[.1,.05,.025])
        for dt in ([.1,.05,.02],[.1,.1,.1],[.1,0,.01],[.1,float('nan'),.025]):
            with self.subTest(dt=dt),self.assertRaises(ValueError):compare(self.fixtures(),dt)

    def test_no_potential_and_partial_presence(self):
        states=self.fixtures()
        for s in states:s.pop('potential_moments')
        result=compare(states,[.1,.05,.025])
        self.assertNotIn('potential_moment',result['variables'])
        self.assertAlmostEqual(result['variables']['gas_energy']['observed_order'],2.,places=9)
        states[0]['potential_moments']=[0.]
        with self.assertRaisesRegex(ValueError,'all states or none'):compare(states,[.1,.05,.025])


if __name__=='__main__':unittest.main()

import unittest
import numpy as np
from predictive_resolution_audit import trace_error


class TraceErrorTests(unittest.TestCase):
    def test_uniform_and_split_reference(self):
        a=np.array([0.,0.]);b=np.array([1.,0.]);normal=np.array([0.,-1.])
        edges=[(0,-1,a,b,normal,1.)]
        left=np.array([[0.,0.],[.5,0.],[.5,1.],[0.,1.]])
        right=left+np.array([.5,0.])
        class Provider:pass
        provider=Provider();provider.last_profiles={(0,0):np.array([1.,3.])}
        result=trace_error(provider,edges,[left,right],[1.,3.])
        self.assertAlmostEqual(result['relative_L1'],0.,places=12)
        provider.last_profiles={(0,0):np.array([2.,2.])}
        result=trace_error(provider,edges,[left,right],[1.,3.])
        self.assertAlmostEqual(result['relative_L1'],.5,places=12)
        square=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        result=trace_error(provider,edges,[square],[2.])
        self.assertAlmostEqual(result['relative_L1'],0.,places=12)


if __name__=='__main__':unittest.main()

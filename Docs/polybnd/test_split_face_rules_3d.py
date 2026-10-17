import unittest
import numpy as np
from split_face_rules_3d import rule,piece_rule


class SharedSplitFaceTests(unittest.TestCase):
    def test_exact_zero_area_piece_and_repeated_cuts(self):
        n=np.array([0.,0,1])
        x,w=piece_rule(np.array([[0.,0,0],[.5,0,0],[1.,0,0]]),n,3)
        self.assertEqual(len(w),0)
        p=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0]])
        x,w=rule(p,n,lambda q:np.column_stack([q[:,0]-.5,.5-q[:,0],q[:,0]-.5]),3)
        self.assertAlmostEqual(w.sum(),1.,places=14)

    def test_common_rule_integrates_all_sample_signs(self):
        p=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0]]);n=np.array([0.,0,1])
        speeds=lambda x:np.column_stack([x[:,0]+x[:,1]-1,x[:,0]-.5,np.zeros(len(x))])
        x,w=rule(p,n,speeds,3)
        self.assertAlmostEqual(w.sum(),1.,places=14)
        self.assertTrue(np.all(w>0))
        np.testing.assert_allclose(w@np.maximum(speeds(x),0),[1/6,1/8,0.],atol=1e-14,rtol=0)
        np.testing.assert_allclose(w@np.minimum(speeds(x),0),[-1/6,-1/8,0.],atol=1e-14,rtol=0)

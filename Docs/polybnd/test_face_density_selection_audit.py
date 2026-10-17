import unittest
from types import SimpleNamespace
import numpy as np
from scipy import sparse
from face_density_selection_audit import select,line_prior_mean,one_sided_interval


class SelectionTests(unittest.TestCase):
    def test_projection_has_known_unique_trace_from_two_starts(self):
        p=dict(k=2,Aub=sparse.csr_matrix((0,4)),bub=np.empty(0),c=np.zeros(4),
            Aeq=sparse.csr_matrix([[1.,1.,0.,0.]]),b=np.ones(1),A=np.ones((1,2)),
            bounds=[(0.,1.),(0.,1.),(0.,0.),(0.,0.)],lower=np.zeros(2),upper=np.ones(2),
            base=SimpleNamespace(x=np.array([1.,0.,0.,0.]),fun=0.),spatial_objective=lambda z:0.)
        for start in (0,1):
            result,z=select(p,np.eye(2),np.array([.25,.75]),np.array([.5,.5]),1e-8,1e-12,start)
            self.assertTrue(result['gate']);np.testing.assert_allclose(z,[.25,.75],atol=1e-10)

    def test_inferred_line_mean_resolves_sliver(self):
        metadata=dict(flat=False,normal=[1.,0.],offset=.993)
        actual=line_prior_mean(np.array([0.,0.]),np.array([1.,0.]),metadata,0.,1.,100.)
        self.assertAlmostEqual(actual,1.+99.*.007,places=12)

    def test_narrow_feasible_range_does_not_stall(self):
        eps=1e-7
        p=dict(k=2,Aub=sparse.csr_matrix((0,4)),bub=np.empty(0),c=np.zeros(4),
            Aeq=sparse.csr_matrix([[1.,1.,0.,0.]]),b=np.ones(1),A=np.ones((1,2)),
            bounds=[(1-eps,1.),(0.,eps),(0.,0.),(0.,0.)],lower=np.array([1-eps,0.]),upper=np.array([1.,eps]),
            base=SimpleNamespace(x=np.array([1.,0.,0.,0.]),fun=0.),spatial_objective=lambda z:0.)
        result,z=select(p,np.eye(2),np.array([.2,.8]),np.array([.5,.5]),1e-8,1e-12,0)
        self.assertTrue(result['gate']);np.testing.assert_allclose(z,[1-eps,eps],atol=1e-10)

    def test_material_interface_has_one_sided_trace(self):
        left=np.array([[0.,0.],[.5,0.],[.5,1.],[0.,1.]])
        right=left+np.array([.5,0.]);a=np.array([.5,0.]);b=np.array([.5,1.]);inward=np.array([-1.,0.])
        for polygon in (left,left[::-1]):self.assertEqual(one_sided_interval(polygon,a,b,inward),(0.,1.))
        self.assertIsNone(one_sided_interval(right,a,b,inward))


if __name__=='__main__':unittest.main()

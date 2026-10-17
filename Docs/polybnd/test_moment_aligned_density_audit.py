import unittest
import numpy as np
from discontinuous_density_audit import clip
from moment_aligned_density_audit import polygon_moments,moment_prior,aligned_parts
from subcell_density_audit import triangle_moments
from face_density_range_audit import trace_weights


class MomentAlignedTests(unittest.TestCase):
    def setUp(self):self.v=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])

    def test_moments_recover_synthetic_oblique_line(self):
        normal=np.array([np.cos(.37),np.sin(.37)]);offset=.7;poly=clip(self.v,normal,offset)
        target=polygon_moments(self.v)+9*polygon_moments(np.asarray(poly.exterior.coords[:-1]))
        neighbor=dict(flat=False,normal=[1.,0.],offset=.5,fraction=.5)
        result=moment_prior(self.v,target,np.array([1.,10.]),neighbor)
        self.assertLess(result['moment_fit_after'],1e-14)
        np.testing.assert_allclose(result['normal'],normal,atol=1e-6)

    def test_aligned_partition_preserves_all_geometry_moments(self):
        meta=dict(flat=False,normal=[1.,0.],offset=.993)
        parts=aligned_parts(self.v,2,meta)
        np.testing.assert_allclose(sum((triangle_moments(v) for v in parts),np.zeros(7)),polygon_moments(self.v),atol=1e-13)
        for part in parts:
            x=part[:,0]-.993
            self.assertFalse(x.min()<-1e-12 and x.max()>1e-12)
        for a,b in zip(self.v,np.roll(self.v,-1,axis=0)):
            weights,_=trace_weights(parts,a,b);self.assertAlmostEqual(weights.sum(),1.,places=10)


if __name__=='__main__':unittest.main()

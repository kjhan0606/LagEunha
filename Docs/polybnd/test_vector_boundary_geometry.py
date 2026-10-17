import unittest
import numpy as np
from shapely.geometry import Polygon
from vector_boundary_geometry import triangle_moments_many,normalized_many,phase_averages
from subcell_density_audit import triangle_moments,normalized
from discontinuous_density_audit import clip
from boundary_layer_feasibility_audit import boundary_parts
from boundary_layer_feasibility_audit import boundary_face_matrix


class VectorGeometryTests(unittest.TestCase):
    def test_triangle_moments(self):
        v=np.random.default_rng(731).uniform(size=(128,3,2));raw=triangle_moments_many(v)
        expected=np.array([triangle_moments(p) for p in v]);np.testing.assert_allclose(raw,expected,rtol=1e-13,atol=1e-15)
        center=np.array([.5,.5]);scale=.2
        np.testing.assert_allclose(normalized_many(raw,center,scale),np.array([normalized(p,center,scale) for p in raw]),rtol=1e-12,atol=1e-14)

    def test_planar_triangle_fractions(self):
        v=np.random.default_rng(42).uniform(size=(128,3,2));normal=np.array([.6,.8]);offset=.55
        result=phase_averages(v,dict(normal=normal,offset=offset),1.,100.)
        expected=np.array([1+99*clip(p,normal,offset).area/Polygon(p).area for p in v])
        np.testing.assert_allclose(result,expected,rtol=1e-10,atol=1e-10)

    def test_thin_layer_integral(self):
        square=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        parts=boundary_parts(square,1e-7,64,moment_function=triangle_moments_many);raw=triangle_moments_many(parts)
        normal=np.array([.6,.8]);offset=.55
        result=phase_averages(parts,dict(normal=normal,offset=offset),1.,100.)
        expected=1+99*clip(square,normal,offset).area
        self.assertAlmostEqual(raw[:,0]@result,expected,places=11)
        self.assertTrue(np.all((result>=1)&(result<=100)))

    def test_aligned_core_preserves_boundary_order(self):
        square=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        prior=dict(flat=False,normal=[.6,.8],offset=.55)
        parts=boundary_parts(square,1e-7,16,moment_function=triangle_moments_many,core_prior=prior)
        count=len(parts)-len(square)*2*16
        self.assertNotEqual(count,4*4**2)
        for a,b in zip(square,np.roll(square,-1,axis=0)):
            w0,w1=boundary_face_matrix(square,len(parts),16,a,b,core_count=count)
            self.assertAlmostEqual(w0.sum(),1.,places=14)
            self.assertAlmostEqual(w1.sum(),0.,places=14)
            self.assertTrue(np.all(w0[:count]==0))


if __name__=='__main__':unittest.main()

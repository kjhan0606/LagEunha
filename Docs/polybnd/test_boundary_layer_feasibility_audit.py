import unittest
import numpy as np
from scipy.optimize import linprog
from subcell_density_audit import triangle_moments
from moment_aligned_density_audit import polygon_moments
from face_moment_feasibility_audit import face_matrix,build,solve
from boundary_layer_feasibility_audit import boundary_parts,trace_interval,boundary_face_matrix


class BoundaryLayerTests(unittest.TestCase):
    def test_tiny_layer_preserves_original_moment_equations(self):
        v=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        parts=boundary_parts(v,1e-7,64);raw=np.array([triangle_moments(p) for p in parts])
        rho=np.array([1.5+p.mean(0)[0] for p in parts]);target=rho@raw
        p=build(v,target,parts,dict(flat=True),target[0],1.,10.);p['moment_row_scale']=1e4
        rows=[]
        for a,b in zip(v,np.roll(v,-1,axis=0)):rows.extend(boundary_face_matrix(v,len(parts),64,a,b))
        F=np.asarray(rows);result=solve(p,F,F@rho/10.)
        self.assertTrue(result['feasible']);self.assertLess(result['moment_error'],1e-10)
        self.assertLess(result['face_error'],1e-9)

    def test_thin_partition_and_trace_moments(self):
        v=np.array([[0.,0.],[1.,0.],[.9,1.],[0.,.8]])
        for radial,segments in ((1e-2,4),(1e-7,64)):
            parts=boundary_parts(v,radial,segments)
            raw=sum((triangle_moments(p) for p in parts),np.zeros(7))
            np.testing.assert_allclose(raw,polygon_moments(v),rtol=1e-10,atol=1e-14)
            for a,b in zip(v,np.roll(v,-1,axis=0)):
                T0,T1=face_matrix(parts,a,b)
                self.assertAlmostEqual(T0.sum(),1.,places=8)
                self.assertAlmostEqual(T1.sum(),0.,places=8)
                self.assertEqual(np.count_nonzero(T0),segments)
                A0,A1=boundary_face_matrix(v,len(parts),segments,a,b)
                np.testing.assert_allclose(A0,T0,atol=1e-8)
                np.testing.assert_allclose(A1,T1,atol=1e-8)
                R0,R1=boundary_face_matrix(v,len(parts),segments,b,a)
                np.testing.assert_allclose(R0,A0);np.testing.assert_allclose(R1,-A1)

    def test_first_moment_interval_agrees_with_lp(self):
        for segments in (4,16,64):
            tau=-1+(2*np.arange(segments)+1)/segments
            for mean in (.02,.5,1.):
                interval=trace_interval(mean,.01,1.01,segments)
                sol=linprog(-tau/segments,A_eq=np.ones((1,segments))/segments,b_eq=[mean],bounds=(.01,1.01),method='highs')
                self.assertTrue(sol.success)
                self.assertAlmostEqual(interval['maximum_first'],-sol.fun,places=12)
                self.assertAlmostEqual(interval['minimum_first'],sol.fun,places=12)


if __name__=='__main__':unittest.main()

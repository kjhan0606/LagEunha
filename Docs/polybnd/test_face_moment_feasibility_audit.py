import unittest
import numpy as np
from subcell_density_audit import triangles,triangle_moments
from face_moment_feasibility_audit import build,solve,face_matrix,reference_face_moments


class FaceFeasibilityTests(unittest.TestCase):
    def setUp(self):self.v=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])

    def test_exact_first_moment_and_reversal(self):
        left=np.array([[0.,0.],[.7,0.],[.7,1.],[0.,1.]])
        right=np.array([[.7,0.],[1.,0.],[1.,1.],[.7,1.]])
        a=np.array([0.,0.]);b=np.array([1.,0.]);normal=np.array([0.,-1.])
        m0,m1,_=reference_face_moments([left,right],[1.,10.],a,b,normal)
        self.assertAlmostEqual(m0,3.7);self.assertAlmostEqual(m1,1.89)
        r0,r1,_=reference_face_moments([left,right],[1.,10.],b,a,normal)
        self.assertAlmostEqual(r0,m0);self.assertAlmostEqual(r1,-m1)

    def test_manufactured_simultaneous_face_constraints(self):
        parts=triangles(self.v,2);rho=np.array([1.5+v.mean(0)[0] for v in parts])
        raw=np.array([triangle_moments(v) for v in parts]);target=rho@raw
        p=build(self.v,target,parts,dict(flat=True),target[0],1.,10.)
        F=[]
        for a,b in zip(self.v,np.roll(self.v,-1,axis=0)):F.extend(face_matrix(parts,a,b))
        F=np.asarray(F);result=solve(p,F,F@rho/10.)
        self.assertTrue(result['feasible']);self.assertLess(result['face_error'],1e-8)
        self.assertLess(result['moment_error'],1e-8)

    def test_impossible_face_bound_is_rejected(self):
        parts=triangles(self.v,2);target=2*sum((triangle_moments(v) for v in parts),np.zeros(7))
        p=build(self.v,target,parts,dict(flat=True),2.,1.,10.)
        F=np.asarray([face_matrix(parts,self.v[0],self.v[1])[0]])
        result=solve(p,F,np.array([2.])) # mean20 exceeds upper10.1
        self.assertFalse(result['feasible']);self.assertEqual(result['classification'],'solver_infeasible')


if __name__=='__main__':unittest.main()

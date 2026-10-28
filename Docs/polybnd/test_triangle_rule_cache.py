import unittest
from unittest.mock import patch
import numpy as np
from power_cells_3d import triangle_rule,face_quadrature


class TriangleRuleCacheTests(unittest.TestCase):
    def test_identical_reference_arithmetic_and_readonly(self):
        for order in (2,4,6,8):
            z,w=np.polynomial.legendre.leggauss(order);z=(z+1)/2;w=w/2
            a,b=np.meshgrid(z,z,indexing='ij');wa,wb=np.meshgrid(w,w,indexing='ij')
            a=a.ravel();b=b.ravel();weights=(wa*wb).ravel()*a
            bary=np.column_stack([1-a,a*(1-b),a*b]);got,q=triangle_rule(order)
            np.testing.assert_array_equal(got,bary);np.testing.assert_array_equal(q,weights)
            self.assertFalse(got.flags.writeable);self.assertFalse(q.flags.writeable)

    def test_face_calls_reuse_only_reference_rule(self):
        triangle_rule.cache_clear();original=np.polynomial.legendre.leggauss
        p=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0]]);n=np.array([0.,0,1])
        with patch('numpy.polynomial.legendre.leggauss',side_effect=original) as gauss:
            x,w=face_quadrature(p,n,4);y,v=face_quadrature(p,n,4)
            self.assertEqual(gauss.call_count,1)
        np.testing.assert_array_equal(x,y);np.testing.assert_array_equal(w,v)
        x[:]=0.;w[:]=0.
        z,q=face_quadrature(p,n,4)
        np.testing.assert_array_equal(z,y);np.testing.assert_array_equal(q,v)


if __name__=='__main__':unittest.main()

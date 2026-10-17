import unittest
import numpy as np
from moment_storage_conditioning import triangle_raw,translated_raw
from moment_aligned_density_audit import polygon_moments


class MomentStorageTests(unittest.TestCase):
    def test_analytic_affine_transform_preserves_raw_moments(self):
        v=np.array([[0.,0.],[1.,0.],[.2,.7]])
        for origin,scale in [(np.array([.45,.6]),.125),(np.array([3.,-2.]),2.)]:
            expected=polygon_moments(v*scale+origin)
            np.testing.assert_allclose(translated_raw(triangle_raw(v),origin,scale),expected,rtol=1e-13,atol=1e-15)

    def test_cap_fraction_is_constructed_in_local_coordinates(self):
        normal=np.array([np.cos(.73),np.sin(.73)])
        for fraction in (1e-2,1e-7,1e-9):
            delta=np.sqrt(2*fraction*normal.prod())
            v=np.array([[1.,1.],[1-delta/normal[0],1.],[1.,1-delta/normal[1]]])
            self.assertLess(abs(float(triangle_raw(v)[0])/fraction-1),1e-10)


if __name__=='__main__':unittest.main()

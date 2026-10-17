import unittest
import numpy as np
from shapely.geometry import Polygon
from shear_interface_audit import moments,exact_on_cells


class InterfaceTests(unittest.TestCase):
    def test_square_moments(self):
        p=Polygon([(0,0),(1,0),(1,1),(0,1)])
        np.testing.assert_allclose(moments(p),[1,.5,.5,1/3,1/4,1/4,1/3],atol=1e-14,rtol=0)

    def test_split_density_reference(self):
        left=np.array([[0,0],[.5,0],[.5,1],[0,1]])
        right=left+np.array([.5,0]);whole=np.array([[0,0],[1,0],[1,1],[0,1]])
        exact,coverage,_=exact_on_cells([whole],[left,right],[1,100],0.)
        self.assertAlmostEqual(exact[0,0],50.5,places=11)
        self.assertLess(coverage,1e-14)

    def test_shear_mass_and_moments(self):
        p=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        t=.1;q=p.copy();q[:,0]+=t*(.2+1.1*q[:,1])
        exact,coverage,shape=exact_on_cells([q],[p],[2.],t)
        self.assertAlmostEqual(exact[0,0],2.,places=12)
        self.assertAlmostEqual(exact[0,1],1.+2*t*(.2+1.1*.5),places=12)
        self.assertLess(coverage+shape,1e-14)

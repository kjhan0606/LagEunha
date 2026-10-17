import unittest
import numpy as np
from split_face_mass_3d import integrate_sample


class SplitFaceMassTests(unittest.TestCase):
    def test_capacity_matches_affine_offset_derivative(self):
        p=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0]]);n=np.array([0.,0,1])
        for sign in (-1,1):
            speed=lambda x:x[:,0]+x[:,1]-.8
            value,capacity=integrate_sample(p,n,speed,lambda x:np.ones(len(x)),sign,4,True)
            h=1e-5
            plus=integrate_sample(p,n,lambda x:speed(x)+h,lambda x:np.ones(len(x)),sign,4)
            minus=integrate_sample(p,n,lambda x:speed(x)-h,lambda x:np.ones(len(x)),sign,4)
            self.assertAlmostEqual((plus-minus)/(2*h),capacity,places=9)

    def test_affine_positive_and_negative_parts(self):
        p=np.array([[0.,0,0],[1,1,0],[1,0,0],[0,1,0]]);n=np.array([0.,0,1])
        speed=lambda x:x[:,0]+x[:,1]-1
        for order in (2,4,6):
            for sign in (-1,1):
                got=integrate_sample(p,n,speed,lambda x:np.ones(len(x)),sign,order)
                self.assertAlmostEqual(got,sign/6,places=14)

    def test_zero_speed_edge_and_empty_halfplane(self):
        p=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0]]);n=np.array([0.,0,1])
        for speed,expected in [(lambda x:x[:,0],.5),(lambda x:x[:,0]*0,0.),(lambda x:x[:,0]-2,0.)]:
            self.assertAlmostEqual(integrate_sample(p,n,speed,lambda x:np.ones(len(x)),1,3),expected,places=14)

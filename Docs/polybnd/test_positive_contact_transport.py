import unittest
import numpy as np
from shapely.geometry import Polygon
from positive_transport_patch import Piece,raw,GAMMA
from positive_contact_transport import capacity,evaluate,compress_materials

class ContactTests(unittest.TestCase):
    def test_density_weighted_capacity(self):
        domain=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        donors=[Piece(Polygon([(0,0),(.5,0),(.5,1),(0,1)]),1.,.9,np.zeros(2),np.zeros((2,2))),Piece(Polygon([(.5,0),(1,0),(1,1),(.5,1)]),4.,.9,np.zeros(2),np.zeros((2,2)))]
        sites=np.array([[.25,.5],[.75,.5]]);target=np.array([1.25,1.25])
        w,cells,it,err=capacity(sites,np.zeros(2),domain,target,donors)
        np.testing.assert_allclose(evaluate(sites,w,domain,donors)[0],target,atol=1e-12)
        self.assertAlmostEqual(cells[0].area,.6875,places=11)

    def test_nonconvex_merge_rejected(self):
        f=lambda p:Piece(Polygon(p),1.,.9,np.zeros(2),np.zeros((2,2)))
        a=f([(0,0),(1,0),(1,1),(0,1)])
        b=f([(2,0),(3,0),(3,1),(2,1)])
        with self.assertRaises(RuntimeError): compress_materials([a,b])

if __name__=='__main__': unittest.main()

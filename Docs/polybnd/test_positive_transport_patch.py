import unittest
import numpy as np
from shapely.geometry import Polygon
from positive_transport_patch import Piece,raw,invariant,move,remap,diagnose,same_field

class PositiveTransportTests(unittest.TestCase):
    def test_free_flight_invariants(self):
        f=Piece(Polygon([(0,0),(1,0),(1,1),(0,1)]),1.4,.9,np.array([.3,-.1]),np.array([[.2,.6],[-.4,-.1]]))
        g=move(f,.03)
        np.testing.assert_allclose(invariant(raw([g])),invariant(raw([f])),atol=3e-15,rtol=3e-15)
        x=np.array([.4,.7]);v=f.b+f.B@x
        np.testing.assert_allclose(g.b+g.B@(x+.03*v),v,atol=3e-16)
        self.assertLess(diagnose([g])['max_R'],1e-28)

    def test_heterogeneous_transfer_keeps_residual(self):
        a=Piece(Polygon([(0,0),(.5,0),(.5,1),(0,1)]),1.,.8,np.array([1.,0.]),np.zeros((2,2)))
        b=Piece(Polygon([(.5,0),(1,0),(1,1),(.5,1)]),2.,1.2,np.array([-1.,.2]),np.zeros((2,2)))
        cell=Polygon([(0,0),(1,0),(1,1),(0,1)])
        groups,audit=remap([a,b],[cell])
        self.assertEqual(len(groups[0]),2)
        np.testing.assert_allclose(raw(groups[0]),raw([a,b]),atol=2e-15,rtol=2e-15)
        self.assertGreater(diagnose(groups[0])['max_R'],.05)
        self.assertGreaterEqual(diagnose(groups[0])['min_R'],-1e-15)
        self.assertFalse(same_field(a,b))

    def test_identical_merge(self):
        a=Piece(Polygon([(0,0),(.5,0),(.5,1),(0,1)]),1.,.8,np.array([.2,0.]),np.array([[0.,1.1],[0.,0.]]))
        b=Piece(Polygon([(.5,0),(1,0),(1,1),(.5,1)]),a.rho,a.u,a.b,a.B)
        groups,audit=remap([a,b],[Polygon([(0,0),(1,0),(1,1),(0,1)])])
        self.assertEqual(len(groups[0]),1)
        self.assertLess(audit['merge_error'],2e-15)

if __name__=='__main__': unittest.main()

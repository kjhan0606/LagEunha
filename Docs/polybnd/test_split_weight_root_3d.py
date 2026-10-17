import unittest
from unittest.mock import patch
from types import SimpleNamespace
import numpy as np
from split_face_mass_3d import evaluate,solve


class SplitWeightRootTests(unittest.TestCase):
    def fixture(self):
        x=np.array([[.25,.5,.5],[.75,.5,.5]])
        face=dict(i=0,j=1,normal=np.array([1.,0,0]),
            vertices=np.array([[.5,0,0],[.5,1,0],[.5,1,1],[.5,0,1]]))
        geometry=SimpleNamespace(build=lambda *a,**k:([None,None],[face],{}))
        raw=np.zeros((2,35));raw[:,0]=1.;L=np.tile(.1*np.eye(3),(2,1,1))
        B=np.array([[0.,.4,.3],[0,0,0],[0,0,0]])
        states=[dict(c=p,v=np.array([v,0.,0.]),B=B,R=.01*np.eye(3),
            density=SimpleNamespace(value=lambda q:np.ones(len(q)))) for p,v in zip(x,[.2,-.1])]
        return geometry,x,raw,L,states

    def test_flux_jacobian_gauge_and_conservation(self):
        g,x,r,L,states=self.fixture();z=np.array([.02,-.03]);direction=np.array([.3,-.2]);h=1e-6
        with patch('split_face_mass_3d.decode',side_effect=states*3):
            flux,J=evaluate(g,x,np.zeros(2),r,L,(np.zeros_like(x),z),return_jacobian=True)
            plus=evaluate(g,x,np.zeros(2),r,L,(np.zeros_like(x),z+h*direction))
            minus=evaluate(g,x,np.zeros(2),r,L,(np.zeros_like(x),z-h*direction))
        np.testing.assert_allclose((plus-minus)/(2*h),J@direction,atol=1e-9,rtol=0)
        np.testing.assert_array_equal(J@np.ones(2),0.)
        self.assertEqual(float(flux.sum()),0.)

    def test_root_solves_without_changing_mass_or_sites(self):
        g,x,r,L,states=self.fixture();original=r.copy();sites=x.copy()
        with patch('split_face_mass_3d.decode',side_effect=lambda p,q:states[int(q[1])]):
            r[1,1]=1.;original=r.copy()
            z,d=solve(g,x,np.zeros(2),r,L,np.zeros_like(x))
        self.assertLess(d['residual'],2e-13);self.assertEqual(z[-1],0.)
        np.testing.assert_array_equal(r,original);np.testing.assert_array_equal(x,sites)

import unittest
import numpy as np
from shapely.geometry import Polygon
from density_closure_audit import rule
from moment_ale_evolve import Density


class ClosureTests(unittest.TestCase):
    def test_rule_moments(self):
        p=Polygon([(0,0),(1,0),(1,1),(0,1)])
        for order in (8,16,32):
            x,w=rule(p,order)
            np.testing.assert_allclose([sum(w),w@x[:,0],w@(x[:,0]*x[:,1])],[1,.5,.25],atol=1e-14,rtol=0)

    def test_constant_fit(self):
        p=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        d=Density(p,2.,np.array([1.,1.]),np.array([[2/3,.5],[.5,2/3]]))
        x,w=rule(Polygon(p),16)
        np.testing.assert_allclose(d.value(x),2.,atol=1e-12,rtol=0)

import unittest
import numpy as np
from numpy.polynomial import Polynomial as P
from quadratic_geometry_gcl import QuadraticPath
from weighted_mass_patch import power_cells
from power_spacetime_geometry import moments
class QuadraticGeometryTests(unittest.TestCase):
    def test_reynolds_moment_derivative(self):
        z=QuadraticPath.__new__(QuadraticPath);z.N=4;z.kind='wave';z.h=0.
        x=np.array([[.24,.23],[.76,.26],[.27,.77],[.73,.74]])
        z.x=[[P([x[i,k],.03*(i-k),.005*(i+1)*(k+1)]) for k in range(2)] for i in range(4)]
        z.w=[P([.003*i,-.002*i,.001*i]) for i in range(4)]
        z.boundary=[[P([a]),P([b]),P([c])] for a,b,c in [(-1,0,0),(1,0,1),(0,-1,0),(0,1,1)]]
        t=.4;h=1e-5;xm,wm=z.state(t-h);xp,wp=z.state(t+h)
        fd=(moments(power_cells(xp,wp,z.domain(t+h)))-moments(power_cells(xm,wm,z.domain(t-h))))/(2*h)
        rate,boundary=z.rate(t)
        np.testing.assert_allclose(rate,fd,atol=3e-10,rtol=0)
        np.testing.assert_allclose(rate.sum(0),boundary,atol=1e-16,rtol=0)
if __name__=='__main__':unittest.main()

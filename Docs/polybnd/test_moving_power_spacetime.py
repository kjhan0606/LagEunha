import unittest
import numpy as np
from moving_power_spacetime import events,integrate
class MovingGeometryTests(unittest.TestCase):
    def test_stationary_limit(self):
        x=np.array([[.25,.25],[.75,.25],[.25,.75],[.75,.75]])
        w=np.array([.02,0.,0.,0.]);dw=-2*w;v=np.zeros_like(x)
        t,_=events(x,v,w,dw);np.testing.assert_allclose(t,[0.,.5,1.],atol=1e-13)
        r=integrate(x,v,w,dw,2)
        self.assertLess(max(r[k] for k in ('area_error','first_moment_error','second_moment_error')),1e-14)
    def test_translating_generators(self):
        x=np.array([[.25,.25],[.75,.25],[.25,.75],[.75,.75]])
        v=np.tile([.03,.01],(4,1));w=np.array([.02,0.,0.,0.]);dw=-2*w
        r=integrate(x,v,w,dw,8)
        self.assertLess(max(r[k] for k in ('area_error','first_moment_error','second_moment_error')),1e-13)
if __name__=='__main__':unittest.main()

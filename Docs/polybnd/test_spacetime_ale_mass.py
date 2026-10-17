import unittest
import numpy as np
from spacetime_ale_mass import D0,velocity,boundary_rows,integrate,endpoint

class MaterialBoundaryTests(unittest.TestCase):
    def test_boundary_equation_and_velocity(self):
        T=.3;rows=boundary_rows(T)
        for t in (.13,.71):
            domain=D0+t*T*velocity(D0)
            for i,row in enumerate(rows):
                p=.37*domain[i]+.63*domain[(i+1)%4];n=np.array([row[0](t),row[1](t)])
                self.assertLess(abs(n@p-row[2](t)),3e-16)
                speed=(row[2].deriv()(t)-np.array([row[0].deriv()(t),row[1].deriv()(t)])@p)/np.linalg.norm(n)
                self.assertAlmostEqual(speed,T*(velocity(p[None,:])[0]@n)/np.linalg.norm(n),places=14)
    def test_shared_flux_identity(self):
        x=np.array([[.25,.25],[.75,.25],[.25,.75],[.75,.75]])
        w=np.array([.01,0.,0.,0.]);dw=np.array([-.02,0.,0.,0.]);T=.03
        flux,b,_=integrate(x,w,dw,T)
        actual=endpoint(x+T*velocity(x),w+dw,T)-endpoint(x,w,0)
        self.assertLess(np.max(abs(actual+flux)),2e-14)
        self.assertLess(np.max(abs(flux.sum(0)-b)),2e-15)
        self.assertLess(abs(b[0]),2e-16)
if __name__=='__main__':unittest.main()

import unittest
import numpy as np
from numpy.polynomial import Polynomial as P
from predicted_spacetime_flux import FluxAudit
from quadratic_geometry_gcl import QuadraticPath
from adaptive_moment_batch import setup
import moment_ale_evolve as base
class FluxTests(unittest.TestCase):
    def test_uniform_rest_mesh_flux_equals_geometric_sweep(self):
        mass,y=setup('wave',2,0.);a=y[:-4].reshape(4,20)
        p=QuadraticPath.__new__(QuadraticPath);p.N=4;p.kind='wave';p.h=.1;p.mass=mass;p.mid_factors=[np.zeros((2,2)) for _ in mass]
        p.x=[[P([v]) for v in row] for row in a[:,:2]];p.w=[P([0.,v]) for v in [0.,.001,-.002,.001]]
        p.boundary=[[P([a]),P([b]),P([c])] for a,b,c in [(-1,0,0),(1,0,1),(0,-1,0),(0,1,1)]]
        audit=FluxAudit.__new__(FluxAudit);audit.path=p;audit.grad=np.zeros((4,2))
        polys,_=base.geometry(a[:,:2],a[:,2],0.);audit.states=[base.decode(poly,np.r_[m,row[3:]]) for poly,m,row in zip(polys,mass,a)]
        f=audit.rate(.3);volume,_=p.rate(.3);flux=f[:-4].reshape(4,4)
        np.testing.assert_allclose(flux[:,0],-volume[:,0],atol=1e-15,rtol=0)
        np.testing.assert_allclose(flux[:,1:3],0.,atol=1e-15,rtol=0)
        np.testing.assert_allclose(flux[:,3],.6/(base.GAMMA-1)*flux[:,0],atol=1e-15,rtol=0)
        np.testing.assert_allclose(flux.sum(0),f[-4:],atol=1e-15,rtol=0)
        derivatives=np.array([w.deriv()(.3) for w in p.w])
        p.w=[P([w(.3)]) for w in p.w]
        np.testing.assert_allclose(audit.rate(.3,derivatives),f,atol=1e-15,rtol=0)
        _,jac=audit.rate(.3,derivatives,return_jacobian=True)
        numerical=np.zeros((4,4));eps=1e-6
        for k in range(4):
            d=np.zeros(4);d[k]=eps
            numerical[:,k]=(audit.rate(.3,derivatives+d)[:-4].reshape(4,4)[:,0]-audit.rate(.3,derivatives-d)[:-4].reshape(4,4)[:,0])/(2*eps)
        np.testing.assert_allclose(jac,numerical,atol=1e-10,rtol=0)
        np.testing.assert_allclose(jac.sum(0),0,atol=1e-15,rtol=0)
if __name__=='__main__':unittest.main()

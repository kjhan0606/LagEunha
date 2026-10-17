import unittest
import numpy as np
from periodic2d_hllc import mesh,initial,rhs,integrate_face,fan,quadrature,exact,GAMMA
from fixedmass_hllc import state,sample


class Periodic2DTests(unittest.TestCase):
    def test_oblique_exact_euler_residual(self):
        points=np.random.default_rng(4).random((20,2)); boost=np.array([.37,-.23]); t=.1; h=1e-5
        def flux(points,axis):
            u=exact(points,t,'wave',boost); rho=u[:,0]; v=u[:,1:3]/rho[:,None]
            p=(GAMMA-1)*(u[:,3]-.5*rho*np.sum(v*v,axis=1))
            f=u*v[:,axis,None]; f[:,axis+1]+=p; f[:,3]+=p*v[:,axis]
            return f
        residual=(exact(points,t+h,'wave',boost)-exact(points,t-h,'wave',boost))/(2*h)
        for axis in range(2):
            shift=np.eye(2)[axis]*h; residual+=(flux(points+shift,axis)-flux(points-shift,axis))/(2*h)
        np.testing.assert_allclose(residual,0.,atol=1e-9)

    def test_area_normal_and_quadrature(self):
        m,y=initial(4,'uniform',np.array([.37,-.23]))
        p,a,e=mesh(y[:,:2],y[:,2],audit=True)
        self.assertAlmostEqual(a.sum(),1.,places=13)
        for poly,area in zip(p,a): self.assertAlmostEqual(quadrature(poly)[1].sum(),area,places=14)
        moved=mesh(y[:,:2]+np.array([.71,-.39]),y[:,2]+.3,audit=True)
        np.testing.assert_allclose(moved[1],a,atol=2e-14)

    def test_area_derivative(self):
        m,y=initial(4,'uniform',np.array([0.,0.]))
        rng=np.random.default_rng(19); v=rng.normal(size=(16,2))*.03; z=rng.normal(size=16)*.001
        x=y[:,:2]; w=y[:,2]+rng.normal(size=16)*.0001
        _,a,edges=mesh(x,w,audit=True); rate=np.zeros(16)
        for i,j,aa,bb,delta,length in edges:
            speed=(delta@v[j]-(v[j]-v[i])@((aa+bb)/2)+(z[i]-z[j])/2)/np.linalg.norm(delta)
            rate[i]+=length*speed; rate[j]-=length*speed
        eps=1e-6
        fd=(mesh(x+eps*v,w+eps*z)[1]-mesh(x-eps*v,w-eps*z)[1])/(2*eps)
        np.testing.assert_allclose(rate,fd,atol=1e-10)

    def test_integrated_flux_and_jacobian(self):
        f=fan(state(1.,[.2,.1],1.),state(.7,[-.1,.3],.8),np.array([1.,0.]))
        flux,density,pieces=integrate_face(f,-2.,2.,.03,.4)
        self.assertEqual(pieces,4)
        eps=1e-6
        derivative=(integrate_face(f,-2.,2.,.03+eps,.4)[0][0]-integrate_face(f,-2.,2.,.03-eps,.4)[0][0])/(2*eps)
        self.assertAlmostEqual(derivative,-density,places=8)
        grid=-2+.03+(np.arange(2000)+.5)*4/2000
        direct=.4*np.mean([sample(f,s)[0] for s in grid],axis=0)
        np.testing.assert_allclose(flux,direct,atol=3e-7)

    def test_translation_and_conservation(self):
        for kind in ('uniform','contact','wave'):
            m,y=initial(4,kind,np.array([.37,-.23])); dy,d=rhs(y,m)
            self.assertLess(d['mass_flux'],2e-13)
            np.testing.assert_allclose(dy[:,3:].sum(0),0.,atol=1e-14)
            np.testing.assert_array_equal(dy[:,:2],y[:,3:5]/m[:,None])
            if kind!='wave':
                np.testing.assert_allclose(dy[:,:2],[ [.37,-.23] ]*16,atol=2e-14)
                np.testing.assert_allclose(dy[:,2:],0.,atol=2e-13)


if __name__=='__main__':unittest.main()

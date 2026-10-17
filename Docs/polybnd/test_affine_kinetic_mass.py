import unittest
import numpy as np
from affine_kinetic_mass import positive_integral,solve,arrays,evaluate,cuts


class ExactMassTests(unittest.TestCase):
    def test_linear_positive_integrals(self):
        a=np.array([1.,-1.,1.,-1.,0.]);b=np.array([3.,-3.,-1.,1.,0.])
        value,fraction=positive_integral(a,b)
        np.testing.assert_allclose(value,[2.,0.,.25,.25,0.]);np.testing.assert_allclose(fraction,[1.,0.,.5,.5,.5])

    def test_weight_root_matches_split_integral(self):
        bins=[dict(i=0,j=1,scale=1.,distance=1.,sl0=np.full(4,-.2),sl1=np.full(4,.8),sr0=np.full(4,-.2),sr1=np.full(4,.8),rl=1.,rr=10.)]
        z,norm,_=solve(bins,2);self.assertLess(norm,2e-13)
        split=cuts(bins[0],z);gx,gw=np.polynomial.legendre.leggauss(2);total=0.
        for lo,hi in zip(split[:-1],split[1:]):
            t=(lo+hi)/2+(hi-lo)*gx/2;s=-.2+t-(z[0]-z[1])/2
            total+=np.dot(gw,np.maximum(s,0)+10*np.minimum(s,0))*(hi-lo)/2
        self.assertAlmostEqual(total,0.,places=12)
        residual,lap=evaluate(arrays(bins),2,z)
        eps=1e-6;plus=evaluate(arrays(bins),2,z+np.array([eps,0.]))[0]
        minus=evaluate(arrays(bins),2,z-np.array([eps,0.]))[0]
        np.testing.assert_allclose((plus-minus)/(2*eps),-lap[:,0],atol=1e-8)

    def test_separating_cold_streams_have_nonunique_weight_root(self):
        bins=[dict(i=0,j=1,scale=1.,distance=1.,sl0=np.full(4,-1.),sl1=np.full(4,-1.),
            sr0=np.full(4,1.),sr1=np.full(4,1.),rl=1.,rr=10.)]
        data=arrays(bins)
        for value in (-1.,0.,1.):
            residual,lap=evaluate(data,2,np.array([value,0.]))
            np.testing.assert_array_equal(residual,0.);np.testing.assert_array_equal(lap,0.)
        self.assertGreater(evaluate(data,2,np.array([-100.,0.]))[0][0],0.)
        self.assertLess(evaluate(data,2,np.array([100.,0.]))[0][0],0.)

    def test_connected_random_graphs(self):
        rng=np.random.default_rng(731)
        for _ in range(20):
            bins=[]
            for i in range(4):
                for j in range(i+1,4):
                    bins.append(dict(i=i,j=j,scale=rng.uniform(.2,1.),distance=rng.uniform(.2,1.),
                        sl0=rng.normal(size=4),sl1=rng.normal(size=4),sr0=rng.normal(size=4),sr1=rng.normal(size=4),
                        rl=rng.uniform(.5,100.),rr=rng.uniform(.5,100.)))
            z,norm,_=solve(bins,4);self.assertLess(norm,2e-13)
            residual,_=evaluate(arrays(bins),4,z)
            self.assertLess(abs(residual.sum()),1e-12)


if __name__=='__main__':unittest.main()

import unittest
import numpy as np
from face_density_range_audit import trace_weights,problem,extremum
from subcell_density_audit import triangles,triangle_moments


class TraceTests(unittest.TestCase):
    def test_face_partition_and_linear_integral(self):
        v=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        for level in (2,3):
            for a,b in zip(v,np.roll(v,-1,axis=0)):
                weights,mid=trace_weights(triangles(v,level),a,b)
                self.assertAlmostEqual(weights.sum(),1.,places=12)
                np.testing.assert_allclose(weights@mid,(a+b)/2,atol=1e-12)

    def test_fixed_uniform_density_has_zero_trace_range(self):
        v=np.array([[0.,0.],[1.,0.],[0.,1.]])
        # With low=high=2, the +/-1% box still permits changes. Set an exact box.
        p=problem(v,2*triangle_moments(v),1,dict(flat=True),2.,2.,2.)
        k=p['k'];z=2*p['raw'][:,0]/p['mass'];p['bounds']=list(zip(z,z))+[(0,None)]*k
        p['lower']=z;p['upper']=z
        weights,_=trace_weights(p['parts'],v[0],v[1]);coefficient=weights*p['mass']/p['raw'][:,0]
        lo=extremum(p,coefficient,1e-8,1.);hi=extremum(p,coefficient,1e-8,-1.)
        self.assertTrue(lo['gate'] and hi['gate'])
        self.assertAlmostEqual(lo['value'],2.,places=10);self.assertAlmostEqual(hi['value'],2.,places=10)


if __name__=='__main__':unittest.main()

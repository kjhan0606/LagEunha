import unittest
from unittest.mock import patch
import numpy as np
from split_face_mass_3d import solve
from test_split_weight_root_3d import SplitWeightRootTests


class ConstantDensity:
    def value(self,points):return np.ones(len(points))


class ParallelSplitMassTests(unittest.TestCase):
    def test_spawn_workers_preserve_root_and_reduction_order(self):
        geo,x,raw,L,states=SplitWeightRootTests().fixture();raw[1,1]=1.
        for state in states:state['density']=ConstantDensity()
        original=geo.build
        def build(*args,**kwargs):
            cells,faces,diag=original(*args,**kwargs)
            return cells,faces*129,diag
        geo.build=build
        with patch('split_face_mass_3d.decode',side_effect=lambda p,r:states[int(r[1])]):
            one,sd=solve(geo,x,np.zeros(2),raw,L,np.zeros_like(x),workers=1)
            two,pd=solve(geo,x,np.zeros(2),raw,L,np.zeros_like(x),workers=2)
        np.testing.assert_array_equal(one,two);self.assertEqual(sd,pd)

    def test_invalid_workers_and_unprepared_parallel_rejected(self):
        geo,x,raw,L,_=SplitWeightRootTests().fixture()
        for workers in (0,-1,1.5):
            with self.subTest(workers=workers),self.assertRaises(ValueError):
                solve(geo,x,np.zeros(2),raw,L,np.zeros_like(x),workers=workers)
        with self.assertRaisesRegex(ValueError,'requires fixed-state'):
            solve(geo,x,np.zeros(2),raw,L,np.zeros_like(x),workers=2,reuse_geometry=False)


if __name__=='__main__':unittest.main()

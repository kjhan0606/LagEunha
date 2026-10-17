import unittest
from unittest.mock import patch
from types import SimpleNamespace
import numpy as np
from quadratic_pressure_reconstruction import prepare,value

class PressureReflectionTests(unittest.TestCase):
    def test_planar_step_has_identical_reconstruction_in_each_row(self):
        n=8;cells=[];states=[]
        for j in range(n):
            for i in range(n):
                cells.append(np.array([[i,j],[i+1,j],[i+1,j+1],[i,j+1]],dtype=float)/n)
                states.append(dict(P=1. if i<n/2 else .1,density=SimpleNamespace(area=1/n**2)))
        for degree in (2,4):
            with patch.dict('os.environ',LAGFORCE_PRESSURE_REFLECT_Y='0,1',LAGFORCE_PRESSURE_DEGREE=str(degree)):
                pp,d=prepare(states,cells)
            self.assertEqual(len(pp),n*n)
            for i in range(n):
                base=value(pp[i],np.array([[i+.2,.2],[i+.8,.8]])/n)
                for j in range(n):
                    v=value(pp[j*n+i],np.array([[i+.2,j+.2],[i+.8,j+.8]])/n)
                    np.testing.assert_allclose(v,base,atol=3e-12,rtol=0)

    def test_invalid_bounds_are_rejected(self):
        with patch.dict('os.environ',LAGFORCE_PRESSURE_REFLECT_Y='1,0'):
            with self.assertRaisesRegex(ValueError,'bounds'):prepare([],[])

if __name__=='__main__':unittest.main()

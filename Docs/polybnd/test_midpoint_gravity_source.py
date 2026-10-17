import unittest
from unittest.mock import patch
import numpy as np
from test_frozen_moment_boundary import state
from gravity_source_2d import apply as original
from midpoint_gravity_source import apply


class MidpointSourceTests(unittest.TestCase):
    def test_affine_acceleration_matches_exact_moment_kick(self):
        m,y,L=state();g=lambda p:.03*p+np.array([.01,-.02]);mask=y[:-4].reshape(-1,20)[:,0]<.3
        with patch.dict('os.environ',LAGFORCE_QUADRATIC_VELOCITY='0'):
            a,A,_=original(0.,y,m,L,g,.001,mask);b,B,d=apply(0.,y,m,L,g,.001,mask)
        np.testing.assert_allclose(a,b,atol=2e-13,rtol=0)
        for x,z in zip(A,B):np.testing.assert_allclose(x@x.T,z@z.T,atol=2e-13,rtol=0)
        np.testing.assert_array_equal(b[:-4].reshape(-1,20)[mask],y[:-4].reshape(-1,20)[mask])

if __name__=='__main__':unittest.main()

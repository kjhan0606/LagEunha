import unittest
from unittest.mock import patch
import numpy as np
from potential_moment_coupled_2d import dense_state


class PotentialMomentPathTests(unittest.TestCase):
    def test_original_potential_endpoints_and_frozen_cells(self):
        r=dict(U=np.array([2.,3.]),du0=np.array([.1,.3]),dum=np.array([.2,.4]),dt=.5,
               frozen_first=np.array([False,True]))
        with patch('gravity_coupled_2d.hydro_dense_state',return_value=('y','L','dy')):
            for s in (0.,.5,1.):
                y,L,dy,U,dU=dense_state(r,np.ones(2),s)
                self.assertEqual(U[1],3.);self.assertEqual(dU[1],0.)
                self.assertAlmostEqual(U[0],2.+.5*(.1*s+.1*s*s))
                self.assertAlmostEqual(dU[0],.1+.2*s)
            self.assertEqual(dense_state(r,np.ones(2),1.)[3][0],2.1)
        np.testing.assert_array_equal(r['U'],[2.,3.])

    def test_reject_extrapolation(self):
        with self.assertRaises(ValueError):dense_state({},[],1.1)

    def test_reject_unknown_gravity_scheme(self):
        from potential_moment_coupled_2d import step
        with self.assertRaisesRegex(ValueError,'unknown gravity'):
            step(0.,None,None,None,None,.001,None,None,None,gravity_scheme='invalid')

if __name__=='__main__':unittest.main()

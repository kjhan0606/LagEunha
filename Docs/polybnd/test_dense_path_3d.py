import unittest
import numpy as np
from moment_state_3d import pack
from moment_evolution_3d import dense_state

class DensePath3DTests(unittest.TestCase):
    def record(self):
        Q=np.full((3,3),.25);np.fill_diagonal(Q,1/3);L=np.eye(3)[None,:,:]*.1
        raw=pack(1.,np.full(3,.5),Q,np.zeros(3),np.zeros((3,3)),L[0]@L[0].T,1.515)[None,:]
        diag=dict(cov_A=[np.zeros((3,3))],cov_K=[np.zeros((3,3))],cov_Qfactor=[np.zeros((3,0))],cov_linearization_factors=L)
        return dict(t=0.,dt=.01,sites=np.full((1,3),.5),weights=np.zeros(1),raw=raw,factors=L,dx0=np.zeros((1,3)),dxmid=np.zeros((1,3)),dw0=np.zeros(1),dwmid=np.zeros(1),dr0=np.zeros((1,35)),drmid=np.zeros((1,35)),diag_mid=diag)

    def test_stationary_primary_factors_and_fixed_mass(self):
        record=self.record();record['dr0'][:,0]=1.;record['drmid'][:,0]=2.
        for s in (0.,.25,.5,1.):
            sites,weights,raw,L,dx,dw=dense_state(record,s)
            np.testing.assert_array_equal(raw,record['raw']);np.testing.assert_array_equal(L,record['factors']);np.testing.assert_array_equal(dx,0.)
        with self.assertRaisesRegex(ValueError,'outside'):dense_state(record,1.1)

    def test_geometry_derivative_is_that_of_recorded_quadratic_path(self):
        record=self.record();record['dx0'][:]=.2;record['dxmid'][:]=-.1;record['dw0'][:]=.3;record['dwmid'][:]=.5
        s=.37;eps=1e-3;center=dense_state(record,s);left=dense_state(record,s-eps);right=dense_state(record,s+eps)
        np.testing.assert_allclose((right[0]-left[0])/(2*eps*record['dt']),center[4],atol=4e-10,rtol=0)
        np.testing.assert_allclose((right[1]-left[1])/(2*eps*record['dt']),center[5],atol=1e-10,rtol=0)

if __name__=='__main__':unittest.main()

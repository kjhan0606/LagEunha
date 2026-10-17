import unittest
from unittest.mock import patch
import numpy as np
from test_frozen_moment_boundary import state
from gravity_coupled_2d import hydro_step,hydro_dense_state,step
from smooth_mass_path import evaluate_smooth_path


class GravityDenseTests(unittest.TestCase):
    def test_record_does_not_change_split_evolution(self):
        m,y,L=state();f=y[:-4].reshape(-1,20)[:,0]<.3
        a,A,d=step(0.,y,m,L,1e-5,lambda p:-.01*p,frozen=f)
        b,B,e=step(0.,y,m,L,1e-5,lambda p:-.01*p,frozen=f,return_record=True)
        np.testing.assert_array_equal(a,b)
        for x,z in zip(A,B):np.testing.assert_array_equal(x,z)
        self.assertIn('hydro_record',e)

    def test_hydro_endpoints_and_prescribed_frozen_path(self):
        m,y,L=state();f=y[:-4].reshape(-1,20)[:,0]<.3
        out,F,d=hydro_step(0.,y,m,L,1e-5,frozen=f,return_record=True)
        record=d['hydro_record']
        for s,expected,expectedL in [(0.,y,L),(1.,out,F)]:
            z,Z,v=hydro_dense_state(record,m,s)
            np.testing.assert_allclose(z,expected,atol=1e-14,rtol=0)
            for a,b in zip(Z,expectedL):np.testing.assert_allclose(a,b,atol=1e-14,rtol=0)
            np.testing.assert_array_equal(z[:-4].reshape(-1,20)[f],y[:-4].reshape(-1,20)[f])
        z,Z,v=hydro_dense_state(record,m,.3);rates=v[:-4].reshape(-1,20)
        with patch('frozen_weight_root.solve',side_effect=AssertionError('audit must not solve root')):
            flux,diag=evaluate_smooth_path(3e-6,z,Z,m,(rates[:,:2],rates[:,2]),frozen_mask=f)
        self.assertTrue(np.all(np.isfinite(flux)))
        bad=rates[:,2].copy();bad[f]=.1
        with self.assertRaisesRegex(ValueError,'nonzero'):
            evaluate_smooth_path(3e-6,z,Z,m,(rates[:,:2],bad),frozen_mask=f)

    def test_mask_transition_is_not_silently_certified(self):
        m,y,L=state();_,_,d=hydro_step(0.,y,m,L,1e-5,return_record=True)
        d['hydro_record']['frozen_middle'][0]=True
        with self.assertRaisesRegex(ValueError,'transition'):
            hydro_dense_state(d['hydro_record'],m,.5)

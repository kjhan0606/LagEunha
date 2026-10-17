import unittest
from unittest.mock import patch
import numpy as np
from test_frozen_moment_boundary import state
from smooth_mass_path import evaluate_smooth_path


class PotentialExchangeTests(unittest.TestCase):
    def test_constant_potential_matches_active_mass_exchange(self):
        m,y,L=state();frozen=y[:-4].reshape(-1,20)[:,0]<.3
        gen=np.zeros((len(m),2));z=np.linspace(0.,.01,len(m));z[frozen]=0.
        with patch.dict('os.environ',LAGFORCE_VECTOR_SMOOTH_PATH='1'):
            old,_=evaluate_smooth_path(0.,y,L,m,(gen,z),frozen_mask=frozen)
            flux,d=evaluate_smooth_path(0.,y,L,m,(gen,z),frozen_mask=frozen,potential=lambda p:np.full(len(p),3.))
        np.testing.assert_allclose(flux,old,atol=1e-14,rtol=0)
        self.assertAlmostEqual(d['active_potential_mass_exchange'],3*flux[~frozen].sum(),places=13)
        self.assertGreater(abs(d['active_potential_mass_exchange']),1e-6)
        np.testing.assert_allclose(d['potential_mass_rates'],3*flux,atol=1e-14,rtol=0)

    def test_internal_faces_cancel_from_total_potential_exchange(self):
        m,y,L=state();gen=np.zeros((len(m),2));z=np.linspace(0.,.01,len(m))
        _,d=evaluate_smooth_path(0.,y,L,m,(gen,z),potential=lambda p:p[:,0]**2+p[:,1])
        self.assertLess(abs(d['active_potential_mass_exchange']),1e-14)
        self.assertLess(abs(d['potential_mass_rates'].sum()),1e-14)

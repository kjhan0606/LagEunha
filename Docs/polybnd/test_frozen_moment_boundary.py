import unittest
from unittest.mock import patch
import numpy as np
from adaptive_moment_batch import setup
from factor_moment_evolve import rebuild
from kinetic_factor_rhs import rhs
from gravity_coupled_2d import step


def state():
    mass,y=setup('wave',3,.01);a=y[:-4].reshape(-1,20);v=np.array([.02,-.01]);a[:,9:11]=mass[:,None]*v
    for row in a:row[11:15]=np.outer(row[3:5],v).ravel()
    L=[np.sqrt(m)*np.diag([.003,.002]) for m in mass];return mass,rebuild(y,mass,L),L


class FrozenBoundaryTests(unittest.TestCase):
    def test_inactive_faces_do_not_call_pressure_solver(self):
        m,y,L=state();mask=np.ones(len(m),dtype=bool)
        with patch.dict('os.environ',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star'):
            with patch('euler_star_pressure.star',side_effect=AssertionError('inactive traction')):
                dy,d=rhs(0.,y,m,covariance_factors=L,frozen_mask=mask)
        np.testing.assert_allclose(dy,0.,atol=1e-14,rtol=0)

    def test_no_frozen_mask_retains_original_rhs(self):
        m,y,L=state();a,d=rhs(0.,y,m,covariance_factors=L);b,e=rhs(0.,y,m,covariance_factors=L,frozen_mask=np.zeros(len(m),dtype=bool))
        np.testing.assert_array_equal(a,b)

    def test_mixed_boundary_preserves_frozen_rates_and_exchange_ledger(self):
        m,y,L=state();frozen=y[:-4].reshape(-1,20)[:,0]<.3
        dy,d=rhs(0.,y,m,covariance_factors=L,frozen_mask=frozen);rates=dy[:-4].reshape(-1,20)
        np.testing.assert_array_equal(rates[frozen],0.)
        total=np.r_[rates[:,9:11].sum(axis=0),rates[:,19].sum(),(rates[:,12]-rates[:,13]).sum()]
        np.testing.assert_allclose(total,dy[-4:],atol=1e-13,rtol=0)
        self.assertLess(d['mass_flux'],2e-13);self.assertEqual(d['generator_velocity_defect'],0.)
        self.assertEqual(d['mass_root_backend'],'frozen_dirichlet_sparse')
        self.assertGreater(np.linalg.norm(d['frozen_exchange_ledger']),0.)

    def test_full_split_step_keeps_frozen_state_and_factors_exact(self):
        m,y,L=state();frozen=y[:-4].reshape(-1,20)[:,0]<.3;calls=[]
        def boundary(t,sites):calls.append(t);return frozen
        out,factors,d=step(0.,y,m,L,1e-5,lambda p:-.01*p,frozen=boundary)
        np.testing.assert_array_equal(out[:-4].reshape(-1,20)[frozen],y[:-4].reshape(-1,20)[frozen])
        for i in np.flatnonzero(frozen):np.testing.assert_array_equal(factors[i],L[i])
        self.assertGreaterEqual(len(calls),4);self.assertLess(d['stage_mass_flux'],2e-13)

    def test_all_frozen_and_invalid_audit_combination(self):
        m,y,L=state();mask=np.ones(len(m),dtype=bool)
        out,factors,d=step(0.,y,m,L,1e-4,lambda p:-p,frozen=mask)
        np.testing.assert_array_equal(out[:-4],y[:-4]);self.assertEqual(d['gravity_work'],0.)
        with self.assertRaisesRegex(ValueError,'explicit'):
            rhs(0.,y,m,covariance_factors=L,frozen_mask=mask,audit_geometry_rates=(np.zeros((len(m),2)),np.zeros(len(m))))

if __name__=='__main__':unittest.main()

import json,unittest
from unittest.mock import patch
import numpy as np
from quadratic_velocity_reconstruction import projected_budget,projected_budget_batch,prepare


class BatchedQuadraticBudgetTests(unittest.TestCase):
    def test_budget_backend_configuration_guard(self):
        from simple_wave_evolution import run_with_executor
        cfg=dict(mass_root_backend='legacy_dense',velocity_reconstruction='affine',quadratic_budget_backend='batched')
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_VELOCITY':'0','LAGFORCE_VECTOR_KINETIC_ROOT':'0','LAGFORCE_VECTOR_QUADRATIC_BUDGET':'0'}),patch('simple_wave_evolution.Path.read_bytes',return_value=json.dumps(cfg).encode()):
            with self.assertRaisesRegex(ValueError,'budget backend configuration'):run_with_executor()

    def test_batch_matches_the_identical_scalar_algorithm_across_rank_and_scale(self):
        rng=np.random.default_rng(417);L=rng.normal(size=(24,2,2));H=rng.normal(size=(24,2,3));G=np.tile(np.eye(3),(24,1,1));mass=np.exp(rng.uniform(-12,2,24))
        L[0]=0.;L[1,:,1]=0.;L[2]*=1e-20;H[2]*=1e-17
        G[:,0,1]=.1;G[:,1,2]=-.2
        for regularized in [False,True]:
            batched=projected_budget_batch(L,mass,H,G,regularized=regularized)
            scalar=[projected_budget(l,m,h,g,regularized=regularized) for l,m,h,g in zip(L,mass,H,G)]
            for j in range(3):np.testing.assert_allclose(batched[j],np.array([v[j] for v in scalar]),atol=3e-12,rtol=3e-11)
            for i in range(24):
                h,w,t=[v[i] for v in batched]
                np.testing.assert_allclose(mass[i]*h@G[i].T@G[i]@h.T+L[i]@t@t.T@L[i].T,L[i]@L[i].T,atol=3e-13,rtol=2e-11)

    def test_whole_reconstruction_batch_matches_scalar_without_state_change(self):
        from test_quadratic_velocity_reconstruction import fixture
        mass,state,L,cells,edges,states=fixture();other=fixture()[-1];before=state.copy()
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_BUDGET':'projected','LAGFORCE_VECTOR_QUADRATIC_BUDGET':'0'}):prepare(states,cells,edges,mass,L)
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_BUDGET':'projected','LAGFORCE_VECTOR_QUADRATIC_BUDGET':'1'}):prepare(other,cells,edges,mass,L)
        for a,b,l,m in zip(states,other,L,mass):
            for k in ('coefficient','standard_coefficient'):np.testing.assert_allclose(a['velocity_quadratic'][k],b['velocity_quadratic'][k],atol=2e-12,rtol=1e-11)
            # The dimensionless residual map contains sqrt(1-sigma**2).
            # Near a fully resolved budget its relative condition number is
            # large; compare the actual physical sample velocity instead.
            np.testing.assert_allclose(l@a['velocity_quadratic']['residual_map']/np.sqrt(m),l@b['velocity_quadratic']['residual_map']/np.sqrt(m),atol=2e-12,rtol=1e-11)
        np.testing.assert_array_equal(state,before)

    def test_full_rhs_and_independent_mass_audit_with_batched_budget(self):
        from test_quadratic_velocity_reconstruction import QuadraticVelocityTests
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_BUDGET':'projected','LAGFORCE_VECTOR_QUADRATIC_BUDGET':'1'}):
            tests=QuadraticVelocityTests();tests.test_mass_only_audit_matches_full_quadratic_rhs();tests.test_positive_factor_equation_matches_raw_covariance_derivative()


if __name__=='__main__':unittest.main()

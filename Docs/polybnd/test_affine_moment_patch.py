import unittest
from affine_moment_patch import affine_case,curved_case


class AffineMomentTests(unittest.TestCase):
    def test_affine_reconstruction_and_ale_moments(self):
        r=affine_case(3)
        for key in ('gradient_error','face_trace_mismatch','generator_point_velocity_error','split_pressure_error','mass_flux_residual','pressure_rate_max','transverse_momentum_rate_max'):
            self.assertLess(r[key],2e-11,key)
        self.assertLess(abs(r['total_energy_rate']),2e-11)
        self.assertLess(abs(r['angular_momentum_rate']),2e-11)
        self.assertLess(min(e['max_raw_moment_rhs_error'] for e in r['finite_difference']),1e-8)
        self.assertGreater(r['naive_pressure_bias'],1e-4)
        self.assertGreater(r['wrong_origin_gradient_error'],1e-3)

    def test_curved_moment_realizability(self):
        r=curved_case(4)
        self.assertGreaterEqual(r['min_moment_covariance_eigenvalue'],-1e-13)
        self.assertGreaterEqual(r['min_actual_ls_residual_covariance_eigenvalue'],-1e-13)
        self.assertLess(r['cross_term_identity_error'],2e-13)
        self.assertLess(r['split_pressure_error'],2e-12)
        self.assertLess(r['spin_error'],2e-12)


if __name__=='__main__':unittest.main()

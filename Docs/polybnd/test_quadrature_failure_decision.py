import unittest
from smooth_path_quadrature import requires_refinement


class QuadratureDecisionTests(unittest.TestCase):
    def setUp(self):
        self.cfg=dict(path_rate_budget=1e-4,generator_path_budget=1e-5,resolve_mass_quadrature_before_reject=True)
        self.checked=dict(ratio=4.,mass_ratio=1.2,physical_ratio=1.2,quadrature_ratio=4.,rate=5e-5,generator=3e-6)

    def test_unresolved_mass_estimate_is_refined_before_a_time_step_rejection(self):
        self.assertTrue(requires_refinement(self.checked,self.cfg))

    def test_resolved_mass_failure_and_pointwise_failures_are_still_rejected(self):
        self.assertFalse(requires_refinement(dict(self.checked,quadrature_ratio=.5),self.cfg))
        self.assertFalse(requires_refinement(dict(self.checked,rate=1.001e-4),self.cfg))
        self.assertFalse(requires_refinement(dict(self.checked,generator=1.001e-5),self.cfg))

    def test_accepted_or_nonadaptive_mesh_is_not_refined(self):
        self.assertFalse(requires_refinement(dict(self.checked,ratio=.9),self.cfg))
        self.assertFalse(requires_refinement(self.checked,self.cfg,False))

    def test_unproven_cost_change_is_opt_in_and_geometry_lift_keeps_its_existing_rule(self):
        cfg=dict(self.cfg,resolve_mass_quadrature_before_reject=False)
        self.assertFalse(requires_refinement(self.checked,cfg))
        self.assertTrue(requires_refinement(self.checked,cfg,integrated_geometry=True))


if __name__=='__main__':unittest.main()

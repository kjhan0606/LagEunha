import unittest
import numpy as np
from density_moment_obstruction import basis,bernstein_matrix,dual_bound,reference_spec,cellwise_supremum_lower_bound


class MomentObstructionTests(unittest.TestCase):
    def setUp(self):self.poly=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])

    def test_bernstein_enclosure_contains_all_sampled_values(self):
        origin=np.array([.5,.5]);scale=np.array([.5,.5]);B=bernstein_matrix(self.poly,origin,scale)
        points=np.random.default_rng(5).uniform(0,1,(1000,2));a=np.array([.2,.4,-.1,.7,.3,-.2])
        self.assertLessEqual(np.max(abs(basis(points,origin,scale)@a)),np.max(abs(B@a))+1e-14)

    def test_constant_density_defect_has_the_exact_mass_lower_bound(self):
        c=.03;delta=np.array([c,0.,0.,c/3,0.,c/3]);value,a=dual_bound(self.poly,delta,[.5,.5],[.5,.5])
        self.assertAlmostEqual(value,c,places=11)

    def test_zero_mass_can_have_a_nonzero_moment_obstruction(self):
        c=.03;delta=np.array([0.,c/3,0.,0.,0.,0.]);value,a=dual_bound(self.poly,delta,[.5,.5],[.5,.5])
        self.assertGreaterEqual(value,c/3-1e-12);self.assertLessEqual(value,c/2+1e-12)

    def test_vortex_density_uses_background_scale_and_own_reference(self):
        function,scale,goal,low,high=reference_spec(dict(problem='steady_vortex',amplitude=.5,config=dict(physical_goals=dict(density_field=.001))))
        from steady_vortex_exact import physical_raw
        self.assertIs(function,physical_raw);self.assertEqual((scale,goal,low,high),(1.,.001,16,24))

    def test_wave_reference_and_unknown_problem_guard(self):
        from simple_wave_exact import physical_raw
        function,scale,goal,low,high=reference_spec(dict(amplitude=.02,config=dict(physical_goal=.01)))
        self.assertIs(function,physical_raw);self.assertEqual((scale,goal,low,high),(.02,.01,16,20))
        with self.assertRaises(ValueError):reference_spec(dict(problem='unknown'))

    def test_local_supremum_bound_uses_physical_area_and_density_scale(self):
        rectangle=self.poly*[2.,1.]+[2000.,-3000.]
        self.assertAlmostEqual(cellwise_supremum_lower_bound(.04,rectangle,.02),1.,places=12)
        with self.assertRaises(ValueError):cellwise_supremum_lower_bound(.04,rectangle,0.)


if __name__=='__main__':unittest.main()

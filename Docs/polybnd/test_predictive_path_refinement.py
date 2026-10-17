import unittest
import numpy as np
from predictive_path_refinement import combine_samples,combine_level
from predictive_evolution_batch import rejected_trial_dt,physical_samples_valid


class RefinedPathTests(unittest.TestCase):
    def test_quadrature_failure_is_separate_from_a_sampled_physical_failure(self):
        r=dict(ratio=100.,rate=1e-11,generator=0.,pressure=.6,factor_R=0.)
        self.assertTrue(physical_samples_valid(r))
        self.assertFalse(physical_samples_valid(dict(r,rate=6e-9)))
        self.assertFalse(physical_samples_valid(dict(r,pressure=-.1)))

    def test_finest_level_integrates_cubic_and_checks_absolute_integral(self):
        def value(s):return dict(mass_rate=[1e-8*(s-.3)**3],generator=0.,pressure=.6,factor_R=0.)
        whole=[value(i/512) for i in range(513)];focus=[value(i/(128*64)) for i in range(129)]
        d=combine_level(whole,focus,np.ones(1),1e-6,1e-4)
        expected=1e-14*((.7**4-.3**4)/4)
        self.assertAlmostEqual(d['signed'][0],expected,places=28)
        self.assertGreaterEqual(d['absolute'][0],d['estimated_absolute_integral'][0])
        self.assertGreaterEqual(d['quad_delta'][0],d['absolute_quadrature_delta'][0])

    def test_focused_replacement_keeps_cubic_integral(self):
        def sample(s):return dict(mass_rate=[1e-8*s**3],generator=0.,pressure=.6,factor_R=0.)
        whole={str(i):sample(i/256) for i in range(257)}
        focus={str(i):sample(i/4096) for i in range(65)}
        d=combine_samples(whole,focus,np.ones(1),1e-6,1e-4)
        self.assertAlmostEqual(d['signed'][0],2.5e-15,places=28)
        self.assertLess(float(d['quad_delta'][0]),1e-28)
        self.assertEqual(d['refined_nodes'],257)

    def test_disabled_search_halves_dt(self):
        dt,state=rejected_trial_dt(.001,1.,{},False)
        self.assertEqual(dt,.0005);self.assertFalse(state.get('used',False))

    def test_alternatives_are_bounded_then_return_to_small_branch(self):
        old=dict(rejections=5);dt,state=rejected_trial_dt(.01,1.,old,True)
        self.assertEqual(dt,.8);self.assertNotIn('used',old)
        dt,state=rejected_trial_dt(dt,1.,state,True);self.assertEqual(dt,.6)
        dt,state=rejected_trial_dt(dt,1.,state,True);self.assertEqual(dt,.4)
        dt,state=rejected_trial_dt(dt,1.,state,True);self.assertEqual(dt,.005)
        dt,state=rejected_trial_dt(dt,1.,state,True);self.assertEqual(dt,.0025)
        self.assertTrue(state['used'])

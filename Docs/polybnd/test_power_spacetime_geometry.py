import unittest
from power_spacetime_geometry import probe

class SpacetimeTests(unittest.TestCase):
    def test_events_and_all_moments(self):
        for r in probe():
            s=r['split']
            self.assertEqual(len(s['events']),3)
            self.assertAlmostEqual(s['events'][1],r['expected_event'],places=12)
            self.assertNotEqual(s['topologies'][0],s['topologies'][1])
            for k in ('area_error','first_moment_error','second_moment_error','global_exchange_error'):
                self.assertLess(s[k],5e-15)
            self.assertGreater(max(r['unsplit'][k] for k in ('area_error','first_moment_error','second_moment_error')),1e-7)
if __name__=='__main__':unittest.main()

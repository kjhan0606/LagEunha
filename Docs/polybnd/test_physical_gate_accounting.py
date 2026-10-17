import unittest
from physical_gate_accounting import retarget_peak,numerical_gates,budget_duration

class GateAccountingTests(unittest.TestCase):
    def fixture(self):
        return dict(absolute=[1e-7],quad_delta=[1e-10]),[1.],dict(rate=2e-5,generator=1e-6,ledger=1e-16,stage_mass=1e-13,min_pressure=.1,min_factor_R=0.),dict(path_mass_budget=1e-5,path_quadrature_budget=1e-7,path_rate_budget=1e-4,generator_path_budget=1e-5)

    def test_longer_duration_does_not_hide_prior_peak(self):
        a,m,p,c=self.fixture();self.assertTrue(all(numerical_gates(a,m,p,c).values()))
        q=retarget_peak(p,.02,.15);self.assertAlmostEqual(q['rate'],1.5e-4);self.assertEqual(p['rate'],2e-5)
        self.assertFalse(numerical_gates(a,m,q,c)['path_rate'])
        self.assertEqual(q['generator'],p['generator'])

    def test_every_historical_state_gate_remains_explicit(self):
        a,m,p,c=self.fixture()
        for key,value in [('rate',.1),('generator',.1),('ledger',1.),('stage_mass',1.),('min_pressure',0.),('min_factor_R',-1e-5)]:
            q=dict(p);q[key]=value;self.assertFalse(all(numerical_gates(a,m,q,c).values()),key)

    def test_invalid_retarget_duration_rejected(self):
        with self.assertRaises(ValueError):retarget_peak({'rate':1.},0.,1.)

    def test_full_target_is_independent_of_intermediate_stop(self):
        self.assertEqual(budget_duration({'end':.2},0.),.2)
        self.assertEqual(budget_duration({'end':.2,'validation_end':10.},0.),10.)
        self.assertEqual(budget_duration({'end':.5,'validation_end':10.},0.),10.)
        self.assertEqual(budget_duration({'end':.5,'validation_end':10.},.1),9.9)

    def test_shorter_or_nonfinite_target_is_rejected(self):
        for cfg in ({'end':.2,'validation_end':.1},{'end':.2,'validation_end':float('nan')},
                    {'end':float('inf')},{'end':0.}):
            with self.assertRaises(ValueError):budget_duration(cfg,0.)

    def test_full_target_tightens_path_mass_and_rate_together(self):
        from smooth_path_quadrature import assess_values
        cfg=dict(end=.2,validation_end=10.,path_mass_budget=1e-5,path_rate_budget=1e-4,
                 path_quadrature_budget=1e-7,generator_path_budget=1e-5)
        nodes=[dict(flux=[1e-6],diag=dict(generator_defect=0.,min_pressure=1.,min_factor_R=0.)) for _ in range(9)]
        short=assess_values(nodes,[1.],.001,.2,cfg)
        full=assess_values(nodes,[1.],.001,budget_duration(cfg,0.),cfg)
        self.assertAlmostEqual(full['rate']/short['rate'],50.)
        self.assertAlmostEqual(full['mass_ratio']/short['mass_ratio'],50.)

if __name__=='__main__':unittest.main()

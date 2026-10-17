import unittest
from physical_gate_accounting import retarget_peak,numerical_gates

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

if __name__=='__main__':unittest.main()

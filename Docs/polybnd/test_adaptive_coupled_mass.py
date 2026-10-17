import unittest
import numpy as np
from adaptive_coupled_mass import simpson,budget_ratio,plain,restore,trial,assess,configuration
from adaptive_moment_batch import setup


class BudgetTests(unittest.TestCase):
    def test_configuration_separates_resolution_and_amplitude(self):
        a=configuration('wave',3,.01,.05,5e-11)
        self.assertNotEqual(a,configuration('wave',4,.01,.05,5e-11))
        self.assertNotEqual(a,configuration('wave',3,.1,.05,5e-11))

    def test_invalid_parameters_rejected(self):
        for n,amp,end in [(2,.01,.05),(3,.6,.05),(3,.01,0.),(3,float('nan'),.05)]:
            with self.assertRaises(ValueError):configuration('wave',n,amp,end,5e-11)

    def test_simpson_cubic_exact(self):
        x=np.linspace(0.,1.,17)
        np.testing.assert_allclose(simpson(np.column_stack([x**3,x**2]),2.),[.5,2/3],atol=1e-15,rtol=0)

    def test_absolute_budget_does_not_cancel(self):
        values=np.linspace(-1.,1.,17)
        self.assertAlmostEqual(float(simpson(values,1.)),0.)
        self.assertGreater(float(simpson(abs(values),1.)),.49)

    def test_budget_additivity(self):
        end=.05;budget=5e-11
        for fraction in (.1,.3,.6):
            self.assertAlmostEqual(budget_ratio(0.,budget*fraction,0.,0.,end*fraction,end),1.)

    def test_serialized_step_reproduces_audit(self):
        mass,y=setup('wave',3,.01);L=[np.zeros((2,2)) for _ in mass]
        record,_,_=trial(0.,y,L,mass,.0001)
        a=assess(record,mass,.05,5);b=assess(restore(plain(record)),mass,.05,5)
        for k in a:np.testing.assert_array_equal(a[k],b[k])


if __name__=='__main__':unittest.main()

"""Static RHS gates only: no time evolution on the login node."""
import unittest
import numpy as np
from fixedmass_ale_evolve import initial, rhs


class RHSTests(unittest.TestCase):
    def test_uniform_rest(self):
        m,y=initial('uniform'); derivative, d=rhs(y,m)
        np.testing.assert_allclose(derivative, 0., atol=1e-13)

    def test_fixed_mass_and_shared_energy_flux(self):
        for kind in ('contact','acoustic'):
            m,y=initial(kind); derivative,d=rhs(y,m)
            self.assertLess(d['mass_flux'],1e-12)
            self.assertLess(abs(derivative[20:24].sum()),1e-12)
            np.testing.assert_allclose(derivative[12:20].reshape(4,2).sum(0),derivative[24:26],atol=1e-13)
            np.testing.assert_allclose(derivative[:8],y[12:20]/np.repeat(m,2),atol=0.)


if __name__ == '__main__':
    unittest.main(verbosity=2)

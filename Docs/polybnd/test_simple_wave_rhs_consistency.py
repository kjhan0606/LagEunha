import unittest
import numpy as np
from simple_wave_evolution import install_domain
from adaptive_moment_batch import setup
from simple_wave_exact import physical_raw
from simple_wave_rhs_consistency import exact_rates
import moment_ale_evolve as model


class SimpleWaveReynoldsTests(unittest.TestCase):
    def test_uniform_field_moving_cell_rates_match_geometric_difference(self):
        mass,state=setup('wave',3,0.);a=state[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],0.)
        gen=np.random.default_rng(12).normal(0.,.01,(9,2));z=np.linspace(-.002,.002,9)
        rate,areas,adot=exact_rates(cells,edges,a[:,:2],gen,z,0.,0.)
        h=1e-6;raw=[]
        for sign in (-1,1):
            polys,_=model.geometry(a[:,:2]+sign*h*gen,a[:,2]+sign*h*z,0.)
            raw.append(np.asarray([physical_raw(p,0.,0.)[0] for p in polys]))
        difference=(raw[1]-raw[0])/(2*h)
        np.testing.assert_allclose(rate,difference,rtol=2e-7,atol=3e-9)
        self.assertLess(abs(rate[:,0].sum()),1e-14)
        np.testing.assert_allclose(adot,rate[:,0],atol=1e-14)


if __name__=='__main__':unittest.main()

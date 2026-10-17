import unittest
import numpy as np
from sod_exact import star,fields,physical_raw,GAMMA


class SodReferenceTests(unittest.TestCase):
    def test_rankine_hugoniot_and_entropy(self):
        s=star();r,u,p,S=s['rhor'],s['u'],s['p'],s['shock']
        left=np.array([r,r*u,p/(GAMMA-1)+.5*r*u*u]);right=np.array([.125,0.,.1/(GAMMA-1)])
        fl=np.array([r*u,r*u*u+p,(left[2]+p)*u]);fr=np.array([0.,.1,0.])
        np.testing.assert_allclose(fl-fr,S*(left-right),rtol=1e-12,atol=1e-13)
        self.assertGreater(p/r**GAMMA,.1/.125**GAMMA)

    def test_split_integrals_conserve_mass_and_energy(self):
        box=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        for t in (0.,.03,.15):
            raw,*_=physical_raw(box,t)
            self.assertAlmostEqual(raw[0],.5625,places=12)
            self.assertAlmostEqual(raw[-1],.55/(GAMMA-1),places=12)
            self.assertAlmostEqual(raw[7],.9*t,places=12)

    def test_wall_guard(self):
        with self.assertRaises(ValueError):fields(np.array([[.1,.5]]),1.)

    def test_initializer_fixed_moments(self):
        import moment_ale_evolve as model
        self.addCleanup(setattr,model,'domain',model.domain)
        self.addCleanup(setattr,model,'boundary_velocity',model.boundary_velocity)
        from sod_evolution import initial,score
        mass,y,L,_=initial(4,1.)
        self.assertAlmostEqual(mass.sum(),.5625,places=13)
        result=score(0.,y,mass,1.,L)
        self.assertLess(result['errors_over_amplitude']['density_field'],1e-11)
        self.assertTrue(result['reference_valid'])

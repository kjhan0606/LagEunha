import unittest
from unittest.mock import patch
import numpy as np
from moment_gravity_kick import moments
from conservative_gravity_kick import kick


class ConservativeGravityTests(unittest.TestCase):
    def test_both_dimensions_preserve_original_mass_position_and_torque(self):
        rng=np.random.default_rng(800)
        for dim in (2,3):
            x=rng.uniform(.5,1.5,(100,dim));dm=rng.uniform(.01,.02,100);v=rng.normal(size=(100,2*dim,dim));raw=moments(x,dm,v,1.)
            new,L,d=kick(raw,x,dm,v,lambda p:-p/np.linalg.norm(p,axis=1)[:,None]**3,.02)
            for k in ('m','q','Q'):np.testing.assert_array_equal(new[k],raw[k])
            np.testing.assert_allclose(new['D']-new['D'].T,raw['D']-raw['D'].T,atol=1e-15,rtol=0)
            self.assertLess(abs(d['thermal_change']),1e-14)
            self.assertGreaterEqual(np.linalg.eigvalsh(L@L.T).min(),0.)
            self.assertAlmostEqual(new['E']-raw['E'],d['gravity_work'],places=14)

    def test_unrepresented_state_is_rejected(self):
        rng=np.random.default_rng(801);x=rng.uniform(.5,1.5,(40,2));dm=np.ones(40)/40;v=np.zeros((40,4,2));raw=moments(x,dm,v,1.);raw['p']=np.array([.1,0.])
        with self.assertRaisesRegex(ValueError,'represent'):kick(raw,x,dm,v,lambda p:p,0.01)

    def test_actual_2d_affine_and_quadratic_state_kicks(self):
        from adaptive_moment_batch import setup
        from factor_moment_evolve import rebuild
        from gravity_source_2d import apply
        from kepler_spline_gravity import acceleration
        for quadratic in ('0','1'):
            with patch.dict('os.environ',LAGFORCE_QUADRATIC_VELOCITY=quadratic,LAGFORCE_QUADRATIC_BUDGET='projected'):
                mass,y=setup('wave',3,.01);L=[np.sqrt(m)*np.diag([.02,.01]) for m in mass];y=rebuild(y,mass,L);old=y[:-4].reshape(-1,20).copy()
                new,factors,d=apply(0.,y,mass,L,lambda p:acceleration(p,center=[0.,0.]),.001);a=new[:-4].reshape(-1,20)
                np.testing.assert_array_equal(a[:,:9],old[:,:9])
                self.assertAlmostEqual(float((a[:,12]-a[:,13]).sum()),float((old[:,12]-old[:,13]).sum()),places=12)
                self.assertLess(d['max_thermal_change'],1e-12)
                self.assertLess(np.max(abs(new-rebuild(new,mass,factors))),1e-13)
                self.assertAlmostEqual(float(a[:,19].sum()-old[:,19].sum()),d['gravity_work'],places=13)

if __name__=='__main__':unittest.main()

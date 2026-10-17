import unittest
import numpy as np
from check_3d_vortex_evolution import fields,initial_raw,velocity_component_integrals
from steady_vortex_exact import fields as fields2d
from moment_evolution_3d import rebuild

class VortexReference3DTests(unittest.TestCase):
    def test_component_integrals_weight_density_and_do_not_cancel(self):
        numerical=np.array([[1.,0.,2.],[-1.,0.,-2.]])
        error,energy=velocity_component_integrals(numerical,np.zeros((2,3)),np.array([2.,4.]),np.array([.25,.75]),.5)
        np.testing.assert_array_equal(error,[2.,0.,4.])
        np.testing.assert_array_equal(energy,[1.75,0.,7.])
        error,energy=velocity_component_integrals(numerical,numerical,np.ones(2),np.array([.25,.75]),1.)
        np.testing.assert_array_equal(error,np.zeros(3))
        with self.assertRaises(ValueError):velocity_component_integrals(numerical,numerical,np.ones(2),np.ones(2),0.)

    def test_extruded_reference_is_stationary_compressible_euler(self):
        x=np.random.default_rng(31).uniform(0,1,(100,3));rho,v,P=fields(x,.3)
        _,_,_,g,dp=fields2d(x[:,:2],amplitude=.3,derivatives=True)
        np.testing.assert_allclose(np.einsum('nij,nj->ni',g,v[:,:2])+dp,0.,atol=2e-14,rtol=0)
        np.testing.assert_allclose(np.trace(g,axis1=1,axis2=2),0.,atol=2e-14,rtol=0)
        np.testing.assert_array_equal(v[:,2],0.);self.assertGreater(P.min(),0.)

    def test_initial_moments_keep_positive_factors_and_thermal_energy(self):
        poly=np.array([[a,b,c] for a in (.3,.7) for b in (.3,.7) for c in (0.,1.)])
        raw,L=initial_raw(poly,.3,14);out=rebuild(raw[None,:],L[None,:,:])[0]
        self.assertAlmostEqual(raw[0],.16,places=13)
        np.testing.assert_allclose(raw,out,atol=3e-14,rtol=0)
        self.assertGreater(out[34]-.5*np.trace(out[25:34].reshape(3,3)),0.)

if __name__=='__main__':unittest.main()

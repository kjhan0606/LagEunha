import unittest
import numpy as np
from adaptive_moment_batch import setup
from kinetic_factor_rhs import rhs
from factor_moment_evolve import rebuild,advance_factors_sylvester as advance
from coupled_factor_path_audit import dense_state


class CoupledTests(unittest.TestCase):
    def test_audit_matches_solved_geometry(self):
        mass,y=setup('wave',3,.01);L=[np.zeros((2,2)) for _ in mass]
        d,diag=rhs(0.,y,mass,covariance_factors=L)
        rates=d[:-4].reshape(len(mass),20)
        audited,ad=rhs(0.,y,mass,covariance_factors=L,
            audit_geometry_rates=(rates[:,:2],rates[:,2]))
        np.testing.assert_array_equal(audited,d)
        self.assertEqual(ad['generator_velocity_defect'],0.)
        self.assertEqual(ad['mass_flux'],diag['mass_flux'])

    def test_dense_endpoints(self):
        mass,y=setup('wave',3,.01);L=[np.zeros((2,2)) for _ in mass]
        y=rebuild(y,mass,L);dt=.001
        d0,d=rhs(0.,y,mass,covariance_factors=L)
        halfL=advance(L,d,dt/2);half=rebuild(y+dt*d0/2,mass,halfL)
        dm,mid=rhs(dt/2,half,mass,covariance_factors=halfL)
        end=rebuild(y+dt*dm,mass,advance(L,mid,dt))
        np.testing.assert_array_equal(dense_state(y,L,d0,dm,mid,mass,dt,0.)[0],y)
        np.testing.assert_allclose(dense_state(y,L,d0,dm,mid,mass,dt,1.)[0],end,rtol=0,atol=1e-15)


if __name__=='__main__':unittest.main()

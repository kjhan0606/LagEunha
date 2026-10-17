import unittest
from unittest.mock import patch
import numpy as np
import moment_ale_evolve as base
from test_frozen_moment_boundary import state
from kinetic_factor_rhs import rhs
from gravity_source_2d import apply
from smooth_mass_path import evaluate_smooth_path


def original(i,poly,raw):return base.Density(poly,raw[0],raw[1:3],raw[3:7].reshape(2,2))


class DensityProviderTests(unittest.TestCase):
    def test_rhs_original_provider_exact(self):
        m,y,L=state();f=y[:-4].reshape(-1,20)[:,0]<.3
        a,d=rhs(0.,y,m,covariance_factors=L,frozen_mask=f)
        b,e=rhs(0.,y,m,covariance_factors=L,frozen_mask=f,density_provider=original)
        np.testing.assert_array_equal(a,b)

    def test_gravity_original_provider_exact(self):
        m,y,L=state();g=lambda p:-.01*p
        a,A,d=apply(0.,y,m,L,g,1e-5)
        b,B,e=apply(0.,y,m,L,g,1e-5,density_provider=original)
        np.testing.assert_array_equal(a,b)
        for x,z in zip(A,B):np.testing.assert_array_equal(x,z)

    def test_path_original_provider_exact(self):
        m,y,L=state();v=np.zeros((len(m),2));z=np.linspace(0.,.001,len(m));phi=lambda p:p[:,0]**3
        a,d=evaluate_smooth_path(0.,y,L,m,(v,z),potential=phi)
        b,e=evaluate_smooth_path(0.,y,L,m,(v,z),potential=phi,density_provider=original)
        np.testing.assert_array_equal(a,b)
        np.testing.assert_array_equal(d['potential_mass_rates'],e['potential_mass_rates'])

if __name__=='__main__':unittest.main()

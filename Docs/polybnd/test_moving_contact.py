import unittest
import numpy as np
from adaptive_moment_batch import setup
from kinetic_factor_rhs import rhs
from moving_contact_check import translated_state


class TranslationTests(unittest.TestCase):
    def test_translation_composes(self):
        mass,y=setup('moving_contact',4,.2);u=np.array([.2,-.1])
        a=translated_state(translated_state(y,mass,u,.02),mass,u,.03)
        np.testing.assert_allclose(a,translated_state(y,mass,u,.05),atol=2e-14,rtol=0)

    def test_initial_rhs(self):
        for speed in (.2,2.):
            mass,y=setup('moving_contact',4,speed);u=np.array([speed,-.5*speed])
            dy,d=rhs(0.,y,mass,covariance_factors=[np.zeros((2,2)) for _ in mass])
            # Centered derivative of a quadratic exact translation is exact up to roundoff.
            h=.01
            exact=(translated_state(y,mass,u,h)-translated_state(y,mass,u,-h))/(2*h)
            np.testing.assert_allclose(dy,exact,atol=1e-10,rtol=0)

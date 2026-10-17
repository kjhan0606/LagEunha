import unittest
import numpy as np
from moment_gravity_kick import moments,kick


class GravityKickTests(unittest.TestCase):
    def test_uniform_acceleration_moments_in_two_and_three_dimensions(self):
        rng=np.random.default_rng(41)
        for dim in (2,3):
            x=rng.normal(size=(100,dim));dm=rng.uniform(.01,.1,100);v=rng.normal(size=(100,2*dim,dim))
            old=moments(x,dm,v,2.);g=np.arange(1,dim+1)*.1;dt=.03
            new,_=kick(x,dm,v,2.,lambda x:np.tile(g,(len(x),1)),dt);boost=dt*g
            np.testing.assert_allclose(new['p'],old['p']+old['m']*boost,atol=1e-14)
            np.testing.assert_allclose(new['D'],old['D']+np.outer(old['q'],boost),atol=1e-14)
            expected=old['Pi']+np.outer(old['p'],boost)+np.outer(boost,old['p'])+old['m']*np.outer(boost,boost)
            np.testing.assert_allclose(new['Pi'],expected,atol=2e-14)
            self.assertAlmostEqual(new['E']-old['E'],old['p']@boost+.5*old['m']*(boost@boost),places=12)

    def test_central_force_preserves_angular_momentum_and_internal_energy(self):
        rng=np.random.default_rng(42)
        for dim in (2,3):
            x=rng.uniform(.5,1.5,(100,dim));dm=np.full(100,.01);v=rng.normal(size=(100,2*dim,dim))
            old=moments(x,dm,v,.1)
            new,_=kick(x,dm,v,.1,lambda x:-x/np.linalg.norm(x,axis=1)[:,None]**3,.02)
            np.testing.assert_allclose(new['D']-new['D'].T,old['D']-old['D'].T,atol=1e-15)
            self.assertAlmostEqual(new['E']-.5*np.trace(new['Pi']),.1,places=14)
            for key in ('m','q','Q'):np.testing.assert_array_equal(new[key],old[key])

    def test_nonlinear_kick_generates_positive_residual_without_clipping(self):
        rng=np.random.default_rng(43);x=rng.uniform(.5,1.5,(150,3));dm=np.full(150,1/150)
        v=np.tile((x@np.diag([.1,.2,.3]))[:,None,:],(1,6,1))
        new,_=kick(x,dm,v,1.,lambda x:-x/np.linalg.norm(x,axis=1)[:,None]**3,.1)
        m,q,Q,p,D,Pi=(new[k] for k in ('m','q','Q','p','D','Pi'))
        S=Q-np.outer(q,q)/m;C=D-np.outer(q,p)/m
        R=Pi-np.outer(p,p)/m-C.T@np.linalg.solve(S,C)
        np.testing.assert_allclose(new['residual_factor']@new['residual_factor'].T,R,atol=1e-15)
        self.assertGreater(np.linalg.eigvalsh(R).min(),0.)

    def test_frozen_position_kick_is_reversible(self):
        rng=np.random.default_rng(44);x=rng.uniform(.5,1.5,(80,2));dm=np.full(80,.02);v=rng.normal(size=(80,4,2))
        gravity=lambda x:-x/np.linalg.norm(x,axis=1)[:,None]**3
        _,forward=kick(x,dm,v,1.,gravity,.1);_,back=kick(x,dm,forward,1.,gravity,-.1)
        np.testing.assert_allclose(back,v,atol=5e-16,rtol=0.)

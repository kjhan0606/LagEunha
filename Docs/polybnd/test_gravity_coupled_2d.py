import unittest
from unittest.mock import patch
import numpy as np
import moment_ale_evolve as base
from adaptive_moment_batch import setup,D0
from gravity_coupled_2d import step,potential_energy


class CoupledGravityTests(unittest.TestCase):
    def test_uniform_freefall_couples_moving_cells_and_potential_energy(self):
        g=np.array([.2,-.1]);end=.02
        with patch.dict('os.environ',LAGFORCE_QUADRATIC_VELOCITY='0',LAGFORCE_QUADRATIC_PRESSURE='0'):
            for dt in (.01,.005):
                mass,y=setup('wave',3,0.);original=y.copy();a0=original[:-4].reshape(-1,20);L=[np.zeros((2,2)) for _ in mass]
                initial_energy=float(a0[:,19].sum())+potential_energy(0.,y,mass,lambda p:-p@g)
                for k in range(round(end/dt)):
                    t=k*dt;half_velocity=(t+dt/2)*g;start_box=D0+.5*t*t*g
                    # Exact boundary of the drift subflow, matching both split endpoints.
                    with patch.object(base,'domain',lambda s:start_box+(s-t)*half_velocity),patch.object(base,'boundary_velocity',lambda p:half_velocity):
                        y,L,d=step(t,y,mass,L,dt,lambda p:np.tile(g,(len(p),1)))
                        self.assertLess(d['stage_mass_flux'],2e-13)
                        energy=float(y[:-4].reshape(-1,20)[:,19].sum())+potential_energy(t+dt,y,mass,lambda p:-p@g)
                        self.assertAlmostEqual(energy,initial_energy,places=12)
                a=y[:-4].reshape(-1,20);shift=.5*end*end*g
                np.testing.assert_allclose(a[:,:2],a0[:,:2]+shift,atol=2e-13,rtol=0)
                np.testing.assert_allclose(a[:,9:11],mass[:,None]*end*g,atol=2e-13,rtol=0)
                np.testing.assert_allclose(a[:,3:5],a0[:,3:5]+mass[:,None]*shift,atol=2e-13,rtol=0)

if __name__=='__main__':unittest.main()

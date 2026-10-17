import unittest
import numpy as np
from compatible_pressure_kick import cloud,faces,kick,fit
from weighted_mass_patch import power_cells

class PressureKickTests(unittest.TestCase):
    def make(self,pressure,velocity):
        sites=np.array([[.22,.24],[.73,.21],[.26,.77],[.76,.74]])
        cells=power_cells(sites,np.zeros(4),np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]]))
        return [cloud(c,1.,pressure,velocity) for c in cells],faces(sites,np.zeros(4),cells)

    def test_constant_pressure_shear_unchanged(self):
        states,e=self.make(lambda x:np.full(len(x),.6),lambda x:np.column_stack([.2+1.1*x[:,1],np.zeros(len(x))]))
        out,d=kick(states,e,.01,boundary='trace')
        for a,b in zip(states,out):
            np.testing.assert_allclose(a.v,b.v,atol=1e-15,rtol=0)
            np.testing.assert_allclose(a.u,b.u,atol=1e-16,rtol=0)
        self.assertLess(d['conservation_error'],1e-15)

    def test_affine_pressure_exact_force(self):
        g=np.array([.2,-.1]);P=lambda x: .8+x@g
        states,e=self.make(P,lambda x:np.zeros_like(x))
        out,d=kick(states,e,.01,boundary='trace',pressure_face=lambda x:.8+x@g)
        for s in out:np.testing.assert_allclose(s.v,np.tile(-.01*g,(len(s.x),1)),atol=2e-16,rtol=0)
        self.assertLess(d['conservation_error'],2e-15)

    def test_compression_pressure_work(self):
        states,e=self.make(lambda x:np.full(len(x),.6),lambda x:-.3*x)
        out,d=kick(states,e,.01,boundary='trace')
        for a,b in zip(states,out):
            self.assertAlmostEqual(b.u.sum()-a.u.sum(),.01*.6*.6*a.polygon.area,places=14)
        self.assertLess(abs(d['kinetic_change']),1e-15)

    def test_closed_boundary_energy_and_finite_residual(self):
        states,e=self.make(lambda x:.6+.2*x[:,0],lambda x:np.column_stack([x[:,1]**2,-.2*x[:,0]**2]))
        before=[s.v.copy() for s in states];out,d=kick(states,e,.003)
        self.assertLess(d['conservation_error'],2e-15)
        self.assertLess(d['residual_covariance_change'],2e-17)
        self.assertGreater(max(np.trace(fit(s)[5]) for s in states),1e-7)
        self.assertEqual(d['ledger'][-1],0.)
        for s,v in zip(states,before): np.testing.assert_array_equal(s.v,v)

    def test_negative_thermal_step_rejected(self):
        states,e=self.make(lambda x:np.full(len(x),.6),lambda x:100*x)
        originals=[s.u.copy() for s in states]
        with self.assertRaisesRegex(RuntimeError,'nonpositive internal energy'):
            kick(states,e,1.,boundary='trace')
        for s,u in zip(states,originals):np.testing.assert_array_equal(s.u,u)

if __name__=='__main__':unittest.main()

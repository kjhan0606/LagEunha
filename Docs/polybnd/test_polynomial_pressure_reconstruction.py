import unittest
from unittest.mock import patch
import numpy as np
from pressure_curvature_kick import fixture,rotating_patch
from periodic2d_hllc import quadrature
from polynomial_pressure_reconstruction import prepare,value,lower_bound,bernstein_transform,powers,monomials
import moment_ale_evolve as model


def pressure_and_gradient(points):
    x,y=(np.asarray(points)-.5).T
    P=.8+.03*x+.02*y+.04*x*x-.01*x*y+.03*y*y+.02*x**3-.015*x*y*y+.01*y**3+.04*x**4+.01*x*x*y*y+.03*y**4
    gx=.03+.08*x-.01*y+.06*x*x-.015*y*y+.16*x**3+.02*x*y*y
    gy=.02-.01*x+.06*y-.03*x*y+.03*y*y+.02*x*x*y+.12*y**3
    return P,np.column_stack([gx,gy])


def data():
    mass,state,L=fixture();a=state[:-4].reshape(-1,20);cells,_=model.geometry(a[:,:2],a[:,2],0.)
    for i,poly in enumerate(cells):
        points,w=quadrature(poly);P,_=pressure_and_gradient(points);a[i,19]=w@P/(model.GAMMA-1)
    states=[model.decode(p,np.r_[m,row[3:]]) for p,m,row in zip(cells,mass,a)]
    return mass,state,L,cells,states


class PolynomialPressureTests(unittest.TestCase):
    def setUp(self):
        # Domain installers mutate module globals; restore them even on failure.
        import moment_ale_evolve as model
        self.addCleanup(setattr, model, 'domain', model.domain)
        self.addCleanup(setattr, model, 'boundary_velocity', model.boundary_velocity)

    def test_degree4_reproduces_the_global_polynomial_and_existing_energy_mean(self):
        mass,state,L,cells,states=data();polys,diag=prepare(states,cells,4);self.assertEqual(diag['pressure_polynomial_limited_cells'],0)
        for i,(s,p,metadata) in enumerate(zip(states,cells,polys)):
            points,w=quadrature(p);P,_=pressure_and_gradient(points)
            np.testing.assert_allclose(value(metadata,points),P,atol=3e-8,rtol=0.)
            self.assertAlmostEqual(float(w@value(metadata,points)),float((model.GAMMA-1)*state[:-4].reshape(-1,20)[i,19]),places=13)

    def test_bernstein_recovery_for_a_degree4_triangle_polynomial(self):
        nodes,inverse,_=bernstein_transform(4);triangle=np.array([[0.,0.],[1.,0.],[0.,1.]]);points=nodes@triangle
        rng=np.random.default_rng(4);metadata=dict(degree=4,center=np.zeros(2),scale=1.,powers=powers(4)[1:],basis_means=np.zeros(14),mean=2.,coefficient=.01*rng.normal(size=14))
        coefficients=inverse@value(metadata,points);p=rng.dirichlet([1.,1.,1.],size=1000)@triangle
        self.assertGreaterEqual(value(metadata,p).min(),coefficients.min()-1e-12);self.assertLessEqual(value(metadata,p).max(),coefficients.max()+1e-12)

    def test_positive_pressure_limiter_preserves_thermal_energy(self):
        mass,state,L,cells,states=data()
        for i,s in enumerate(states):s['P']=1e-5 if i==len(states)//2 else .5+(i%3)
        polys,diag=prepare(states,cells,4);self.assertGreater(diag['pressure_polynomial_limited_cells'],0)
        for s,p,metadata in zip(states,cells,polys):
            self.assertGreaterEqual(lower_bound(metadata,p),.01*metadata['mean']-2e-11)
            points,w=quadrature(p);self.assertAlmostEqual(float(w@value(metadata,points)),s['P']*s['density'].area,places=12)

    def test_pressure_force_and_position_momentum_for_the_quartic_patch(self):
        from kinetic_factor_rhs import rhs
        mass,state,L,cells,states=data()
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_PRESSURE':'1','LAGFORCE_PRESSURE_DEGREE':'4','LAGFORCE_QUADRATIC_VELOCITY':'0','LAGFORCE_VECTOR_KINETIC_ROOT':'1','LAGFORCE_KINETIC_SMOOTHING':'1'}):dy,diag=rhs(0.,state,mass,covariance_factors=L)
        da=dy[:-4].reshape(-1,20)
        for i,(poly,s) in enumerate(zip(cells,states)):
            if np.any(s['c']<.2) or np.any(s['c']>.8):continue
            points,w=quadrature(poly);_,g=pressure_and_gradient(points)
            np.testing.assert_allclose(da[i,9:11]/mass[i],-(w@g)/mass[i],atol=3e-7,rtol=0.)
            np.testing.assert_allclose(da[i,11:15].reshape(2,2)/mass[i],-(points.T@(w[:,None]*g))/mass[i],atol=3e-7,rtol=0.)
        self.assertEqual(diag['pressure_reconstruction'],'polynomial_cell_mean_degree4');self.assertLess(diag['mass_flux'],1e-12);self.assertEqual(diag['generator_velocity_defect'],0.)

    def test_degree4_keeps_the_rotating_euler_patch_moment_identity(self):
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_PRESSURE':'1','LAGFORCE_PRESSURE_DEGREE':'4','LAGFORCE_QUADRATIC_VELOCITY':'0','LAGFORCE_VECTOR_KINETIC_ROOT':'1','LAGFORCE_KINETIC_SMOOTHING':'1'}):r=rotating_patch()
        self.assertLess(r['maximum_raw_moment_rate_error_per_mass'],3e-8)


if __name__=='__main__':unittest.main()

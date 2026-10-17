import json,unittest
from unittest.mock import patch
import numpy as np
from pressure_curvature_kick import fixture,inspect,HESSIAN
from periodic2d_hllc import quadrature
from quadratic_pressure_reconstruction import prepare,value,polygon_statistics,lower_bernstein_bound,bernstein_lower_certificate
import moment_ale_evolve as model


class QuadraticPressureTests(unittest.TestCase):
    def setUp(self):
        # Domain installers mutate module globals; restore them even on failure.
        import moment_ale_evolve as model
        self.addCleanup(setattr, model, 'domain', model.domain)
        self.addCleanup(setattr, model, 'boundary_velocity', model.boundary_velocity)

    def test_cheap_bound_encloses_detailed_bernstein_coefficients(self):
        states,cells=self.states();polys,_=prepare(states,cells);rng=np.random.default_rng(19)
        for metadata,p in zip(polys,cells):
            metadata=dict(metadata,coefficient=rng.normal(size=5))
            self.assertLessEqual(bernstein_lower_certificate(metadata,p),lower_bernstein_bound(metadata,p)+1e-14)

    def test_cheap_certificate_does_not_change_any_pressure_coefficients(self):
        states,cells=self.states()
        for i,s in enumerate(states):s['P']=1e-5 if i==len(states)//2 else 1.+i%3
        fast,fastdiag=prepare(states,cells)
        with patch('quadratic_pressure_reconstruction.bernstein_lower_certificate',return_value=-float('inf')):reference,diag=prepare(states,cells)
        self.assertEqual(fastdiag,diag)
        for a,b in zip(fast,reference):np.testing.assert_array_equal(a['coefficient'],b['coefficient']);self.assertEqual(a['positivity_scale'],b['positivity_scale'])

    def states(self):
        mass,y,L=fixture();a=y[:-4].reshape(-1,20);cells,_=model.geometry(a[:,:2],a[:,2],0.)
        return [model.decode(p,np.r_[m,row[3:]]) for p,m,row in zip(cells,mass,a)],cells

    def test_polygon_statistics_are_translation_stable(self):
        poly=np.array([[0.,0.],[2.,0.],[2.,1.],[0.,1.]])
        area,center,cov=polygon_statistics(poly);a,c,s=polygon_statistics(poly+[2000.,-3000.])
        self.assertEqual(area,2.);np.testing.assert_allclose(center,[1.,.5],atol=1e-15)
        np.testing.assert_allclose(cov,np.diag([1/3,1/12]),atol=1e-15)
        np.testing.assert_allclose(s,cov,atol=1e-15);np.testing.assert_allclose(c,center+[2000.,-3000.],atol=1e-12);self.assertEqual(a,area)

    def test_global_quadratic_pressure_and_exact_geometric_cell_means(self):
        states,cells=self.states();polynomials,diag=prepare(states,cells)
        self.assertEqual(diag['pressure_quadratic_limited_cells'],0)
        for s,p,metadata in zip(states,cells,polynomials):
            points,w=quadrature(p);z=points-.5;exact=.8+.5*np.einsum('ni,ij,nj->n',z,HESSIAN,z)
            np.testing.assert_allclose(value(metadata,points),exact,atol=2e-10,rtol=0.)
            self.assertAlmostEqual(float(w@value(metadata,points)/w.sum()),s['P'],places=13)

    def test_bernstein_limiter_preserves_mean_and_positive_pressure(self):
        states,cells=self.states()
        for i,s in enumerate(states):s['P']=1e-5 if i==len(states)//2 else 1.+i%3
        polys,diag=prepare(states,cells);self.assertGreater(diag['pressure_quadratic_limited_cells'],0)
        for s,p,metadata in zip(states,cells,polys):
            self.assertGreaterEqual(lower_bernstein_bound(metadata,p),.01*s['P']-2e-14)
            points,w=quadrature(p);self.assertAlmostEqual(float(w@value(metadata,points)/w.sum()),s['P'],places=12)

    def test_quadratic_pressure_produces_the_exact_affine_acceleration_gradient(self):
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_PRESSURE':'1','LAGFORCE_QUADRATIC_VELOCITY':'0','LAGFORCE_VECTOR_KINETIC_ROOT':'1','LAGFORCE_KINETIC_SMOOTHING':'1'}):r=inspect()
        self.assertEqual(r['pressure_reconstruction'],'quadratic_cell_mean')
        self.assertLess(r['maximum_mean_acceleration_error'],1e-9);self.assertLess(r['maximum_gradient_error'],3e-8)
        self.assertLess(r['mass_residual'],1e-12);self.assertEqual(r['generator_defect'],0.)

    def test_quadratic_pressure_shared_flux_conserves_global_ledgers(self):
        from test_quadratic_velocity_reconstruction import fixture as velocity_fixture
        from kinetic_factor_rhs import rhs
        mass,state,L,*_=velocity_fixture()
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_PRESSURE':'1','LAGFORCE_QUADRATIC_VELOCITY':'1','LAGFORCE_QUADRATIC_BUDGET':'projected','LAGFORCE_VECTOR_KINETIC_ROOT':'1'}):dy,diag=rhs(0.,state,mass,covariance_factors=L)
        a=dy[:-4].reshape(-1,20)
        np.testing.assert_allclose(a[:,9:11].sum(0),dy[-4:-2],atol=3e-12,rtol=0.)
        self.assertAlmostEqual(float(a[:,19].sum()),dy[-2],places=12)
        self.assertAlmostEqual(float(np.sum(a[:,12]-a[:,13])),dy[-1],places=12)
        self.assertLess(diag['mass_flux'],1e-12);self.assertEqual(diag['generator_velocity_defect'],0.)

    def test_rotating_euler_patch_has_correct_all_moment_reynolds_rates(self):
        from pressure_curvature_kick import rotating_patch
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_PRESSURE':'1','LAGFORCE_QUADRATIC_VELOCITY':'0','LAGFORCE_VECTOR_KINETIC_ROOT':'1','LAGFORCE_KINETIC_SMOOTHING':'1'}):r=rotating_patch()
        self.assertGreater(r['interior_cells'],0);self.assertLess(r['maximum_raw_moment_rate_error_per_mass'],3e-10)
        self.assertLess(r['mass_residual'],1e-12);self.assertEqual(r['generator_defect'],0.)

    def test_config_rejects_a_different_runtime_pressure_closure(self):
        from simple_wave_evolution import run_with_executor
        config=dict(mass_root_backend='legacy_dense',velocity_reconstruction='affine',pressure_reconstruction='quadratic_cell_mean')
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_PRESSURE':'0','LAGFORCE_QUADRATIC_VELOCITY':'0','LAGFORCE_VECTOR_KINETIC_ROOT':'0'}),patch('simple_wave_evolution.Path.read_bytes',return_value=json.dumps(config).encode()):
            with self.assertRaisesRegex(ValueError,'pressure reconstruction configuration'):run_with_executor()


if __name__=='__main__':unittest.main()

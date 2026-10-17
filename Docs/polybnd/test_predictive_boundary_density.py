import unittest
import numpy as np
from adaptive_moment_batch import setup
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic
import reference_kinetic_factor_rhs as reference
from predictive_boundary_density import BoundaryDensity


class PredictiveDensityTests(unittest.TestCase):
    def test_initial_contact_profile_preserves_numerical_density(self):
        mass,y=setup('contact',3,0.);a=y[:-4].reshape(-1,20);polys,edges=model.geometry(a[:,:2],a[:,2],0.)
        raw=np.column_stack([mass,a[:,3:]]);provider=BoundaryDensity();trace,diag=provider(polys,edges,raw)
        self.assertLess(diag['closure_moment_error'],1e-8)
        self.assertGreater(diag['closure_min_density'],0.)
        for edge_id,(i,j,*_) in enumerate(edges):
            for cell in (i,j):
                if cell<0:continue
                mean=mass[cell]/model.Density(polys[cell],mass[cell],a[cell,3:5],a[cell,5:9].reshape(2,2)).area
                for t in (.001,.5,.999):self.assertAlmostEqual(trace(cell,edge_id,t),mean,places=9)

    def test_mass_path_matches_full_split_rhs(self):
        from predictive_mass_path import evaluate_path
        mass,y=setup('shear',3,0.);L=[np.zeros((2,2)) for _ in mass];a=y[:-4].reshape(-1,20)
        velocity=np.column_stack([.2+1.1*a[:,1],np.zeros(len(mass))]);z=.001*np.sin(np.arange(len(mass)))
        provider=BoundaryDensity(segments=4)
        _,diag=kinetic.rhs(0.,y,mass,covariance_factors=L,audit_geometry_rates=(velocity,z),
            face_density_provider=provider,face_segments=4,face_order=2,exact_bin_transport=True)
        residual,fast=evaluate_path(0.,y,L,mass,(velocity,z),provider)
        np.testing.assert_allclose(residual,diag['mass_rates'],atol=1e-11)

        vector=BoundaryDensity(segments=4,vector_geometry=True)
        residual_vector,_=evaluate_path(0.,y,L,mass,(velocity,z),vector)
        np.testing.assert_allclose(residual_vector,diag['mass_rates'],atol=1e-11)
        self.assertLess(fast['generator_defect'],1e-12)

    def test_default_kinetic_rhs_is_preserved(self):
        mass,y=setup('shear',3,0.);L=[np.zeros((2,2)) for _ in mass]
        expected,de=reference.rhs(0.,y,mass,covariance_factors=L)
        actual,da=kinetic.rhs(0.,y,mass,covariance_factors=L)
        np.testing.assert_allclose(actual,expected,rtol=0,atol=1e-14)
        self.assertAlmostEqual(da['mass_flux'],de['mass_flux'],places=14)

    def test_mass_path_with_nonzero_covariance_matches_full_rhs(self):
        from predictive_mass_path import evaluate_path
        from factor_moment_evolve import rebuild
        mass,y=setup('shear',3,0.)
        L=[np.array([[.003,.001],[-.0002,.002]])*np.sqrt(m) for m in mass]
        y=rebuild(y,mass,L);a=y[:-4].reshape(-1,20)
        velocity=np.column_stack([.2+1.1*a[:,1],np.zeros(len(mass))]);z=.001*np.cos(np.arange(len(mass)))
        provider=BoundaryDensity(segments=4)
        _,diag=kinetic.rhs(0.,y,mass,covariance_factors=L,audit_geometry_rates=(velocity,z),
            face_density_provider=provider,face_segments=4,face_order=4,exact_bin_transport=True)
        residual,_=evaluate_path(0.,y,L,mass,(velocity,z),provider)
        np.testing.assert_allclose(residual,diag['mass_rates'],atol=1e-11)
        vector=BoundaryDensity(segments=4,vector_geometry=True)
        residual_vector,_=evaluate_path(0.,y,L,mass,(velocity,z),vector)
        np.testing.assert_allclose(residual_vector,diag['mass_rates'],atol=1e-11)


if __name__=='__main__':unittest.main()

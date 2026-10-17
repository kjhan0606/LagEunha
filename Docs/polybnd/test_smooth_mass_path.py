import unittest
import tempfile,json
from pathlib import Path
import numpy as np
from adaptive_moment_batch import setup
from kinetic_factor_rhs import rhs
from smooth_mass_path import evaluate_smooth_path,sample_record,sample_record_file,initialize_worker


class SmoothMassPathTests(unittest.TestCase):
    def test_continuous_density_mass_matches_full_rhs_audit(self):
        mass,y=setup('wave',3,.001);L=[np.zeros((2,2)) for _ in mass]
        dy,diag=rhs(0.,y,mass,covariance_factors=L)
        rates=dy[:-4].reshape(-1,20);gen=rates[:,:2]+np.array([.0002,-.0001])
        z=rates[:,2]+np.linspace(-.001,.001,len(mass))
        _,full=rhs(0.,y,mass,covariance_factors=L,audit_geometry_rates=(gen,z))
        flux,audit=evaluate_smooth_path(0.,y,L,mass,(gen,z))
        np.testing.assert_allclose(flux,full['mass_rates'],rtol=1e-12,atol=1e-14)
        self.assertAlmostEqual(audit['generator_defect'],.0002,places=14)
        self.assertAlmostEqual(audit['min_pressure'],full['min_pressure'],places=14)

    def test_record_worker_hash_and_same_sample(self):
        from adaptive_coupled_mass import atomic
        from predictive_evolution_batch import record_digest
        mass,y=setup('wave',3,.001);zeros=np.zeros_like(y);L=[np.zeros((2,2)) for _ in mass]
        record=dict(t=0.,dt=.001,y=y,L=L,d0=zeros,dm=zeros,
            diag=dict(cov_A=L,cov_K=L,cov_linearization_factors=L,cov_Qfactor=L))
        direct=sample_record(record,mass,.375)
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'pending.json';atomic(p,dict(record=record))
            initialize_worker();cached=sample_record_file(str(p),record_digest(record),mass,.375)
            np.testing.assert_array_equal(direct['flux'],cached['flux']);self.assertEqual(direct['diag'],cached['diag'])
            with self.assertRaises(ValueError):sample_record_file(str(p),'0'*64,mass,.375)

    def test_shared_pressure_gradient_reproduces_a_linear_field(self):
        from types import SimpleNamespace
        from kinetic_factor_rhs import pressure_gradients
        centers=np.array([[.2,.2],[.8,.2],[.2,.8],[.8,.8]])
        states=[dict(density=SimpleNamespace(origin=c),P=.6+.2*c[0]-.3*c[1]) for c in centers]
        edges=[(0,1),(0,2),(1,3),(2,3)]
        origin,grad=pressure_gradients(states,edges)
        np.testing.assert_array_equal(origin,centers);np.testing.assert_allclose(grad,np.tile([.2,-.3],(4,1)),atol=1e-14)


if __name__=='__main__':unittest.main()

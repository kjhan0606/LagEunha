import unittest
from unittest.mock import patch
import numpy as np
from kinetic_speed_smoothing import split
from vector_kinetic_root import record_arrays,evaluate
from kinetic_factor_rhs import evaluate_records,rhs


class KineticSpeedSmoothingTests(unittest.TestCase):
    def test_positive_split_identity_and_derivatives(self):
        s=np.array([-1e6,-.2,0.,.2,1e6]);p,m,dp,dm=split(s,1e-3)
        self.assertTrue(np.all(p>0));self.assertTrue(np.all(m<0))
        np.testing.assert_allclose(p+m,s,atol=1e-14,rtol=0.)
        np.testing.assert_allclose(dp+dm,1.,atol=3e-16,rtol=0.)
        for x in [-.03,0.,.03]:
            p,m,dp,dm=split(np.array([x]),.005);h=1e-7
            a,b,*_=split(np.array([x+h]),.005);c,d,*_=split(np.array([x-h]),.005)
            np.testing.assert_allclose((a-c)/(2*h),dp,atol=1e-9,rtol=0.)
            np.testing.assert_allclose((b-d)/(2*h),dm,atol=1e-9,rtol=0.)

    def test_smooth_graph_jacobian_matches_records_and_finite_difference(self):
        rng=np.random.default_rng(83);records=[]
        for i,j in [(0,1),(1,2),(2,0)]:
            records.append((i,j,1.,.01,.3,rng.normal(.1,.3,(4,2)),rng.normal(.1,.3,(4,2)),1.,1.,np.array([1.,0.]),np.zeros(2),1.,1.,1.,0.,.03))
        data=record_arrays(records);z=np.array([.02,-.01,0.]);r,J,_=evaluate(data,3,z);old,L,_=evaluate_records(records,3,z)
        np.testing.assert_allclose(r,old,atol=1e-14,rtol=0.);np.testing.assert_allclose(J.toarray(),L,atol=1e-14,rtol=0.)
        self.assertGreater(np.linalg.eigvalsh(J.toarray()[:-1,:-1]).min(),0.)
        for k in range(3):
            d=np.eye(3)[k]*1e-7;p,*_=evaluate(data,3,z+d,False);m,*_=evaluate(data,3,z-d,False)
            np.testing.assert_allclose(-(p-m)/2e-7,J.toarray()[:,k],atol=2e-9,rtol=0.)

    def test_full_positive_covariance_identity_and_mass_audit(self):
        from test_quadratic_velocity_reconstruction import QuadraticVelocityTests,fixture
        from geometry_mass_lift import rates
        with patch.dict('os.environ',{'LAGFORCE_KINETIC_SMOOTHING':'1','LAGFORCE_QUADRATIC_BUDGET':'projected','LAGFORCE_QUADRATIC_VELOCITY':'1','LAGFORCE_VECTOR_KINETIC_ROOT':'1'}):
            QuadraticVelocityTests().test_mass_only_audit_matches_full_quadratic_rhs()
            QuadraticVelocityTests().test_positive_factor_equation_matches_raw_covariance_derivative()
            mass,state,L,*_=fixture();dy,diag=rhs(0.,state,mass,covariance_factors=L);g,gd=rates(0.,state,L,mass)
        np.testing.assert_allclose(g,dy[:-4].reshape(-1,20)[:,:3],atol=2e-12,rtol=0.)
        self.assertEqual(diag['kinetic_smoothing_coefficient'],1.)

    def test_affine_constant_density_patch_is_preserved(self):
        from adaptive_moment_batch import setup
        from factor_moment_evolve import rebuild
        mass,state=setup('shear',3,0.);L=[np.zeros((2,2)) for _ in mass];state=rebuild(state,mass,L)
        with patch.dict('os.environ',{'LAGFORCE_KINETIC_SMOOTHING':'0','LAGFORCE_QUADRATIC_VELOCITY':'0','LAGFORCE_VECTOR_KINETIC_ROOT':'1'}):old,_=rhs(0.,state,mass,covariance_factors=L)
        with patch.dict('os.environ',{'LAGFORCE_KINETIC_SMOOTHING':'1','LAGFORCE_QUADRATIC_VELOCITY':'0','LAGFORCE_VECTOR_KINETIC_ROOT':'1'}):new,diag=rhs(0.,state,mass,covariance_factors=L)
        np.testing.assert_allclose(new,old,atol=2e-11,rtol=0.);self.assertLess(diag['mass_flux'],2e-13)

    def test_incompatible_saved_model_is_rejected(self):
        from simple_wave_evolution import run_with_executor
        with patch('pathlib.Path.read_bytes',return_value=b'{"end":0.1}'),patch.dict('os.environ',{'LAGFORCE_KINETIC_SMOOTHING':'1','LAGFORCE_VECTOR_KINETIC_ROOT':'0','LAGFORCE_QUADRATIC_VELOCITY':'0'}):
            with self.assertRaisesRegex(ValueError,'kinetic smoothing configuration'):run_with_executor()


if __name__=='__main__':unittest.main()

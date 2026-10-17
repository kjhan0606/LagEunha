import unittest
import numpy as np
from physical_field_accuracy_audit import errors,merge,finish_norms,environment,reconstructed_raw


class PhysicalFieldAuditTests(unittest.TestCase):
    def test_high_quadrature_raw_energy_counts_residual_velocity_variance_once(self):
        points=np.array([[0.,0.],[1.,0.],[0.,1.],[1.,1.]]);w=np.full(4,.25);rho=np.full(4,2.);v=np.tile([.3,-.2],(4,1));P=np.full(4,.6);R=np.array([[.01,.003],[.003,.02]])
        r=reconstructed_raw(points,w,rho,v,P,R,5/3);np.testing.assert_allclose(r[:7],[2.,1.,1.,1.,.5,.5,1.],atol=1e-15);np.testing.assert_allclose(r[7:9],[.6,-.4],atol=1e-15);np.testing.assert_allclose(r[13:17].reshape(2,2),2*(np.outer(v[0],v[0])+R),atol=1e-15);self.assertAlmostEqual(r[17],1.06,places=14)

    def test_reconstructed_momentum_and_energy_transform_with_a_velocity_boost(self):
        rng=np.random.default_rng(97);points=rng.random((20,2));w=rng.random(20);rho=1+rng.random(20);v=rng.normal(size=(20,2));P=1+rng.random(20);R=np.diag([.03,.07]);boost=np.array([.4,-.6]);r=reconstructed_raw(points,w,rho,v,P,R,5/3);b=reconstructed_raw(points,w,rho,v+boost,P,R,5/3)
        np.testing.assert_allclose(b[7:9],r[7:9]+r[0]*boost,atol=2e-14);np.testing.assert_allclose(b[9:13].reshape(2,2),r[9:13].reshape(2,2)+np.outer(r[1:3],boost),atol=2e-14);self.assertAlmostEqual(b[17],r[17]+boost@r[7:9]+.5*r[0]*np.dot(boost,boost),places=12)

    def test_optional_numerical_backends_remain_explicit_in_field_audits(self):
        cfg=dict(power_cell_evaluation_backend='seeded_active',mass_path_evaluation_backend='cell_batched',moment_flux_evaluation_backend='batched');e=environment(cfg)
        for k in ['LAGFORCE_ACTIVE_POWER_CELLS','LAGFORCE_SEEDED_POWER_CELLS','LAGFORCE_VECTOR_SMOOTH_PATH','LAGFORCE_VECTOR_MOMENT_FLUX']:self.assertEqual(e[k],'1')
        self.assertEqual(environment(dict(power_cell_evaluation_backend='exhaustive'))['LAGFORCE_ACTIVE_POWER_CELLS'],'0')
        with self.assertRaises(ValueError):environment(dict(power_cell_evaluation_backend='unknown'))

    def test_integrated_scalar_vector_norms_and_scales(self):
        left=dict(density=np.array([1.,3.]),velocity=np.array([[3.,4.],[0.,0.]]));right={k:np.zeros_like(v) for k,v in left.items()}
        rows=finish_norms(errors(left,right,[.5,.5],dict(density=2.,velocity=.5)))
        self.assertEqual(rows['density']['L1'],1.);self.assertAlmostEqual(rows['density']['L2'],np.sqrt(1.25));self.assertEqual(rows['velocity']['maximum_sampled'],10.)

    def test_merging_preserves_integrals_and_maximum(self):
        rows={};merge(rows,dict(density=dict(L1=.1,L2_squared=.02,maximum_sampled=.3)));merge(rows,dict(density=dict(L1=.2,L2_squared=.03,maximum_sampled=.2)))
        self.assertAlmostEqual(rows['density']['L1'],.3);self.assertAlmostEqual(rows['density']['L2_squared'],.05);self.assertEqual(rows['density']['maximum_sampled'],.3)

    def test_model_environment_is_explicit_and_unknown_closure_rejected(self):
        e=environment(dict(velocity_reconstruction='quadratic_budget_projected',pressure_reconstruction='quadratic_cell_mean',quadratic_budget_backend='batched'))
        self.assertEqual(e['LAGFORCE_QUADRATIC_PRESSURE'],'1');self.assertEqual(e['LAGFORCE_VECTOR_QUADRATIC_BUDGET'],'1');self.assertEqual(environment({})['LAGFORCE_QUADRATIC_VELOCITY'],'0')
        with self.assertRaises(ValueError):environment(dict(velocity_reconstruction='unknown'))
        e=environment(dict(pressure_reconstruction='polynomial_cell_mean_degree4'))
        self.assertEqual(e['LAGFORCE_PRESSURE_DEGREE'],'4');self.assertEqual(e['LAGFORCE_QUADRATIC_PRESSURE'],'1')
        with self.assertRaises(ValueError):environment(dict(pressure_reconstruction='unknown'))


if __name__=='__main__':unittest.main()

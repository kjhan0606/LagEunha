import unittest
import numpy as np
from common_interface_contact import initialize
from adaptive_moment_batch import setup
from predictive_boundary_density import BoundaryDensity
import moment_ale_evolve as model


class GlobalMomentDirectionTests(unittest.TestCase):
    def test_common_direction_is_inferred_and_local_moments_preserved(self):
        mass,y,_,_=initialize(3,100.,23.);a=y[:-4].reshape(-1,20)
        cells,edges=model.geometry(a[:,:2],a[:,2],0.);raw=np.column_stack([mass,a[:,3:]]);before=raw.copy()
        p=BoundaryDensity('global_moment',vector_geometry=True,align_core=True,construct_prior_witness=True,
            smooth_phase_floor=True,phase_fraction_floor=1e-7,refine_moment_direction=True)
        _,d=p(cells,edges,raw)
        direction=np.array([np.cos(np.radians(23.)),np.sin(np.radians(23.))])
        self.assertLess(float(np.linalg.norm(np.asarray(p.last_global_prior['normal'])-direction)),1e-10)
        self.assertLess(d['closure_moment_error'],1e-8);self.assertGreater(d['closure_min_density'],0.)
        np.testing.assert_array_equal(before,raw)

    def test_uniform_state_needs_no_interface(self):
        mass,y=setup('shear',3,0.);a=y[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],0.)
        p=BoundaryDensity('global_moment',vector_geometry=True,align_core=True,construct_prior_witness=True,
            smooth_phase_floor=True,phase_fraction_floor=1e-7,refine_moment_direction=True)
        _,d=p(cells,edges,np.column_stack([mass,a[:,3:]]))
        self.assertLess(d['closure_moment_error'],1e-8)
        self.assertTrue(all(x['flat'] for x in p.last_priors))

    def test_split_face_order2_matches_order4_for_affine_patch(self):
        import kinetic_factor_rhs as kinetic
        mass,y,_,_=initialize(3,100.,23.);L=[np.zeros((2,2)) for _ in mass]
        p=BoundaryDensity('global_moment',vector_geometry=True,align_core=True,construct_prior_witness=True,
            smooth_phase_floor=True,phase_fraction_floor=1e-7,refine_moment_direction=True)
        args=dict(covariance_factors=L,face_density_provider=p,face_segments=64,exact_bin_transport=True)
        dy2,_=kinetic.rhs(0.,y,mass,face_order=2,**args)
        dy4,_=kinetic.rhs(0.,y,mass,face_order=4,**args)
        np.testing.assert_allclose(dy2,dy4,rtol=0,atol=1e-11)

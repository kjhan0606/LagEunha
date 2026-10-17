import unittest
import numpy as np
from neighbor_density_audit import phase_floor_blend,interface_prior
from face_density_selection_audit import line_prior_mean
from face_moment_feasibility_audit import build
from subcell_density_audit import triangles
from moment_aligned_density_audit import polygon_moments


class PhaseFloorTests(unittest.TestCase):
    def test_finite_resolution_floor_is_explicit_and_dimensionless(self):
        from predictive_boundary_density import BoundaryDensity
        low=1.;high=100.;epsilon=1e-7;threshold=epsilon*(high-low)
        self.assertEqual(phase_floor_blend(low+.5*threshold,low,high,epsilon),0.)
        self.assertAlmostEqual(phase_floor_blend(low+2.5*threshold,low,high,epsilon),.5,places=10)
        self.assertEqual(phase_floor_blend(low+4.1*threshold,low,high,epsilon),1.)
        with self.assertRaises(ValueError):BoundaryDensity(phase_fraction_floor=epsilon)
        with self.assertRaises(ValueError):BoundaryDensity(smooth_phase_floor=True,phase_fraction_floor=-1.)

    def test_blend_endpoints_and_small_endpoint_slopes(self):
        low=1.;high=100.;threshold=32*np.finfo(float).eps*high
        self.assertEqual(phase_floor_blend(low+.5*threshold,low,high),0.)
        self.assertEqual(phase_floor_blend(low+threshold,low,high),0.)
        self.assertAlmostEqual(phase_floor_blend(low+2.5*threshold,low,high),.5,places=12)
        self.assertEqual(phase_floor_blend(low+4*threshold,low,high),1.)
        self.assertLess(phase_floor_blend(low+(1+1/64)*threshold,low,high),1e-4)
        self.assertGreater(phase_floor_blend(low+(4-1/64)*threshold,low,high),1-1e-4)

    def test_blended_face_and_volume_objective(self):
        square=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        metadata=dict(flat=False,normal=[1.,0.],offset=.5,blend=.5)
        value=line_prior_mean(np.array([1.,0.]),np.array([1.,1.]),metadata,6.,1.,11.)
        self.assertEqual(value,8.5)
        p=build(square,6*polygon_moments(square),triangles(square,1),metadata,6.,1.,11.)
        self.assertAlmostEqual(p['spatial_objective'](np.full(p['k'],6/11)),2.5/6,places=12)

    def test_prior_is_unchanged_away_from_floor(self):
        square=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        centers=np.array([[0.,.5],[.5,.5],[1.,.5]]);means=np.array([1.,50.,100.]);neighbors=[{1},{0,2},{1}]
        original,p0=interface_prior(1,square,centers,means,neighbors)
        smooth,p1=interface_prior(1,square,centers,means,neighbors,smooth_phase_floor=True)
        self.assertEqual(p0,p1)
        np.testing.assert_array_equal(original(square),smooth(square))
        finite,p2=interface_prior(1,square,centers,means,neighbors,smooth_phase_floor=True,phase_fraction_floor=1e-7)
        self.assertEqual(p0,p2)
        np.testing.assert_array_equal(original(square),finite(square))

    def test_finite_blend_keeps_moment_direction_and_blend_metadata(self):
        from unittest.mock import patch
        from common_interface_contact import initialize
        import predictive_boundary_density as density
        import moment_ale_evolve as model
        mass,y,_,_=initialize(3,100.,23.);a=y[:-4].reshape(-1,20)
        cells,edges=model.geometry(a[:,:2],a[:,2],0.);raw=np.column_stack([mass,a[:,3:]])
        original=density.interface_prior
        def blended(*args,**kwargs):
            value,metadata=original(*args,**kwargs)
            return value,metadata if metadata['flat'] else dict(metadata,blend=.5)
        provider=density.BoundaryDensity('moment',vector_geometry=True,align_core=True,
            construct_prior_witness=True,smooth_phase_floor=True,phase_fraction_floor=1e-7)
        with patch.object(density,'interface_prior',side_effect=blended),patch.object(density,'moment_prior',wraps=density.moment_prior) as fitted:
            _,diag=provider(cells,edges,raw)
        self.assertGreater(fitted.call_count,0)
        self.assertEqual(fitted.call_count,sum(not p['flat'] for p in provider.last_priors))
        for p in provider.last_priors:
            if not p['flat']:
                self.assertEqual(p['blend'],.5)
                self.assertIn('moment_fit_after',p)
        self.assertLess(diag['closure_moment_error'],1e-8)


if __name__=='__main__':unittest.main()

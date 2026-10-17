import unittest
import numpy as np
from adaptive_moment_batch import setup
from moment_aligned_density_audit import polygon_moments
from curved_contact_rhs import affine_state_from_position_moments
from shear_contact_check import affine_errors
import moment_ale_evolve as model


class CurvedAffineStateTests(unittest.TestCase):
    def test_saved_rhs_loader_accepts_rhs_payload_not_dense_record(self):
        from curved_contact_first_step import load_saved_rhs
        data=dict(dy=[0.,1.],diag=dict(mass_flux=0.,cov_linearization_factors=[[[0.,0.],[0.,0.]]],
            cov_A=[[[0.,1.],[0.,0.]]],cov_K=[[[0.,0.],[0.,0.]]],cov_Qfactor=[[[],[]]]))
        restored=load_saved_rhs(data)
        self.assertEqual(restored['diag']['cov_Qfactor'][0].shape,(2,0))
        self.assertEqual(restored['diag']['mass_flux'],0.)
        np.testing.assert_array_equal(restored['dy'],[0.,1.])

    def test_position_moments_and_affine_velocity_are_preserved(self):
        mass,skeleton=setup('shear_contact',3,100.);a=skeleton[:-4].reshape(-1,20)
        cells,_=model.geometry(a[:,:2],a[:,2],0.)
        geometry=np.asarray([polygon_moments(v) for v in cells])
        rho=np.linspace(1.,100.,len(mass));raw=geometry*rho[:,None]
        m,y,L=affine_state_from_position_moments(skeleton,raw,geometry)
        np.testing.assert_array_equal(m,raw[:,0]);b=y[:-4].reshape(-1,20)
        np.testing.assert_array_equal(b[:,3:9],raw[:,1:7])
        np.testing.assert_array_equal(b[:,:3],a[:,:3])
        self.assertLess(max(affine_errors(y,m).values()),1e-12)
        for poly,mi,row in zip(cells,m,b):
            state=model.decode(poly,np.r_[mi,row[3:]])
            self.assertAlmostEqual(state['P'],.6,places=12)
            self.assertLess(np.max(abs(state['R'])),1e-12)


if __name__=='__main__':unittest.main()

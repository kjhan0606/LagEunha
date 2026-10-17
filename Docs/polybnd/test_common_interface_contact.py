import unittest
import numpy as np
from common_interface_contact import initialize
from shear_contact_check import affine_errors


class CommonInitializationTests(unittest.TestCase):
    def test_common_mass_and_affine_moments(self):
        for n in (3,4,8):
            for angle in (0.,23.):
                mass,y,_,_=initialize(n,100.,angle)
                self.assertAlmostEqual(float(mass.sum()),50.5,places=11)
                self.assertLess(max(affine_errors(y,mass).values()),1e-12)

    def test_direct_prior_witness_preserves_moments_and_budget(self):
        import moment_ale_evolve as model
        from predictive_boundary_density import BoundaryDensity
        from vector_boundary_geometry import triangle_moments_many
        mass,y,_,_=initialize(3,100.,23.);a=y[:-4].reshape(-1,20)
        polys,edges=model.geometry(a[:,:2],a[:,2],0.);raw=np.column_stack([mass,a[:,3:]])
        provider=BoundaryDensity('moment',vector_geometry=True,align_core=True,construct_prior_witness=True)
        _,d=provider(polys,edges,raw)
        self.assertEqual(d['closure_prior_witnesses'],len(mass));self.assertEqual(d['closure_lp_calls'],0)
        for i,((parts,rho),fit) in enumerate(zip(provider.last_cells,provider.last_fits)):
            reconstructed=rho@triangle_moments_many(parts)
            self.assertLess(float(np.max(abs(reconstructed-raw[i,:7]))/mass[i]),1e-10)
            self.assertLessEqual(fit['objective'],1e-5)
            self.assertGreater(float(rho.min()),0.)


if __name__=='__main__':unittest.main()

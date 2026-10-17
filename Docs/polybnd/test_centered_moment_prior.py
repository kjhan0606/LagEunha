import os,unittest
from unittest.mock import patch
import numpy as np
from discontinuous_density_audit import clip
from moment_aligned_density_audit import polygon_moments,moment_prior
from centered_moment_prior import fit_centered_prior


class CenteredMomentPriorTests(unittest.TestCase):
    def test_planar_moments_in_shifted_scaled_cells(self):
        for shift,scale in ((np.zeros(2),1.),(np.array([.45,.6]),.08)):
            vertices=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])*scale+shift
            for angle in (0.,.4):
                normal=np.array([np.cos(angle),np.sin(angle)])
                offset=float(normal@(shift+np.array([.65,.5])*scale))
                cut=clip(vertices,normal,offset)
                target=polygon_moments(vertices)+99*polygon_moments(np.asarray(cut.exterior.coords[:-1]))
                seed=dict(flat=False,normal=[1.,0.])
                result=fit_centered_prior(vertices,target,np.array([1.,100.]),seed,refine=True)
                self.assertLess(np.linalg.norm(np.asarray(result['normal'])-normal),2e-10)
                self.assertLess(abs(result['offset']-offset),2e-11)

    def test_explicit_dispatch_does_not_reuse_legacy_cache(self):
        vertices=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        target=polygon_moments(vertices)+99*polygon_moments(np.array([[.5,0.],[1.,0.],[1.,1.],[.5,1.]]))
        args=(vertices,target,np.array([1.,100.]),dict(flat=False,normal=[1.,0.]))
        with patch.dict(os.environ,{'LAGFORCE_CENTERED_MOMENT_FIT':'0'}):old=moment_prior(*args)
        with patch.dict(os.environ,{'LAGFORCE_CENTERED_MOMENT_FIT':'1'}):new=moment_prior(*args)
        self.assertNotIn('centered_moment_fit',old)
        self.assertTrue(new['centered_moment_fit'])

    def test_uniform_phase_stays_flat(self):
        v=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        self.assertTrue(fit_centered_prior(v,polygon_moments(v),np.array([1.,100.]),dict(flat=True))['flat'])


if __name__=='__main__':unittest.main()

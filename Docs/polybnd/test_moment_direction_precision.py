import unittest
import numpy as np
from discontinuous_density_audit import clip
from moment_aligned_density_audit import polygon_moments,_fit_moment_prior


class MomentDirectionPrecisionTests(unittest.TestCase):
    def test_polish_reduces_residual_without_using_reference_direction(self):
        vertices=np.array([[.1,.2],[.9,.16],[.95,.9],[.05,.85]])
        for angle in (0.,23.):
            n=np.array([np.cos(np.radians(angle)),np.sin(np.radians(angle))]);offset=n@np.array([.5,.5])
            high=clip(vertices,n,offset);target=polygon_moments(vertices)+99*polygon_moments(np.asarray(high.exterior.coords[:-1]))
            means=np.array([1.,100.]);seed=dict(flat=False,normal=[1.,0.],offset=.5)
            old=_fit_moment_prior(vertices,target,means,seed)
            new=_fit_moment_prior(vertices,target,means,seed,refine=True)
            self.assertLessEqual(new['moment_fit_after'],old['moment_fit_after']+1e-30)
            self.assertLess(float(np.linalg.norm(np.asarray(new['normal'])-n)),1e-10)
            self.assertIn('moment_polish_iterations',new)

    def test_default_is_the_unpolished_fit(self):
        vertices=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        target=polygon_moments(vertices)+99*polygon_moments(np.array([[.5,0.],[1.,0.],[1.,1.],[.5,1.]]))
        args=(vertices,target,np.array([1.,100.]),dict(flat=False,normal=[.8,.6],offset=.5))
        self.assertEqual(_fit_moment_prior(*args),_fit_moment_prior(*args,refine=False))


if __name__=='__main__':unittest.main()

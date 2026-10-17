import unittest
import numpy as np
from area_prior_density_audit import interface_averages
from subcell_density_audit import triangles,triangle_moments,fit


class AreaPriorTests(unittest.TestCase):
    def test_thin_corner_missed_by_centroid(self):
        v=np.array([[0.,0.],[1.,0.],[0.,1.]])
        seed=lambda x:np.where(x[:,0]>=.99,100.,1.)
        meta=dict(flat=False,normal=[1.,0.],offset=.99)
        self.assertEqual(seed(np.array([v.mean(0)]))[0],1.)
        self.assertAlmostEqual(interface_averages([v],seed,meta,1.,100.)[0],1.0099,places=12)

    def test_partition_prior_mass(self):
        v=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        seed=lambda x:np.where(x[:,0]>=.993,100.,1.)
        meta=dict(flat=False,normal=[1.,0.],offset=.993)
        for level in (1,2,3):
            parts=triangles(v,level);areas=np.array([triangle_moments(p)[0] for p in parts])
            means=interface_averages(parts,seed,meta,1.,100.)
            self.assertAlmostEqual(means@areas,1.+99.*.007,places=12)

    def test_fit_uses_averages_callback(self):
        v=np.array([[0.,0.],[1.,0.],[0.,1.]])
        result,_,rho=fit(v,2*triangle_moments(v),1,lambda x:np.ones(len(x)),
            density_bounds=(2.,2.),part_averages=lambda parts:np.full(len(parts),2.))
        self.assertTrue(result['moment_gate'] and result['bound_gate'])
        self.assertLess(result['prior_mass_relative_error'],1e-12)
        np.testing.assert_allclose(rho,2.,atol=1e-12)


if __name__=='__main__':unittest.main()

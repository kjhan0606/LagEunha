import unittest
from unittest.mock import patch
import numpy as np
from kinetic_factor_rhs import evaluate_records,solve,rhs
from vector_kinetic_root import record_arrays,evaluate,solve_records
from simple_wave_evolution import initial


class VectorKineticRootTests(unittest.TestCase):
    def records(self):
        rng=np.random.default_rng(82);records=[]
        for i,j in [(0,1),(1,2),(2,3),(3,0)]:
            normal=np.array([1.,0.]);vl=rng.normal(.1,.3,(4,2));vr=rng.normal(.1,.3,(4,2))
            records.append((i,j,1.,.01,.3,vl,vr,1.+i/3,1.+j/3,normal))
        return records

    def test_vector_residual_jacobian_and_fluxes_match_records(self):
        records=self.records();z=np.array([.01,-.02,.03,0.])
        a,L,values=evaluate_records(records,4,z);b,J,(wl,wr,wn)=evaluate(record_arrays(records),4,z)
        np.testing.assert_allclose(a,b,rtol=1e-14,atol=1e-14);np.testing.assert_allclose(L,J.toarray(),rtol=1e-14,atol=1e-14)
        for index,(left,right,speed) in enumerate(values):
            np.testing.assert_allclose(wl[index],left,rtol=1e-14,atol=1e-14);np.testing.assert_allclose(wr[index],right,rtol=1e-14,atol=1e-14)
            self.assertAlmostEqual(wn[index],speed,places=15)

    def test_same_monotone_mass_root(self):
        records=self.records()
        with patch.dict('os.environ',{'LAGFORCE_VECTOR_KINETIC_ROOT':'0'}):old,_,rold=solve(records,4)
        new,_,rnew=solve_records(records,4)
        self.assertLess(rold,2e-13);self.assertLess(rnew,2e-13)
        np.testing.assert_allclose(new,old,rtol=1e-11,atol=1e-12)

    def test_separate_gauges_use_minimum_norm_step(self):
        records=[];n=np.array([1.,0.]);v=np.tile([.2,0.],(4,1))
        for i,j in [(0,1),(2,3)]:records.append((i,j,1.,0.,1.,v,v,1.,1.,n))
        with patch.dict('os.environ',{'LAGFORCE_VECTOR_KINETIC_ROOT':'0'}):old,_,_=solve(records,4)
        new,_,residual=solve_records(records,4)
        np.testing.assert_allclose(new,old,rtol=1e-12,atol=1e-13);self.assertLess(residual,2e-13)

    def test_full_nonlinear_wave_rhs_changes_only_roundoff(self):
        mass,state,L,_=initial(3,.02)
        with patch.dict('os.environ',{'LAGFORCE_VECTOR_KINETIC_ROOT':'0'}):old,do=rhs(0.,state,mass,covariance_factors=L)
        with patch.dict('os.environ',{'LAGFORCE_VECTOR_KINETIC_ROOT':'1'}):new,dn=rhs(0.,state,mass,covariance_factors=L)
        np.testing.assert_allclose(new,old,rtol=1e-10,atol=1e-11)
        self.assertLess(dn['mass_flux'],2e-13);self.assertEqual(dn['mass_root_backend'],'vector_sparse')
        for a,b in zip(do['cov_Qfactor'],dn['cov_Qfactor']):np.testing.assert_allclose(a@a.T,b@b.T,rtol=1e-9,atol=1e-12)


if __name__=='__main__':unittest.main()

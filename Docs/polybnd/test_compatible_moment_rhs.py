import unittest
import numpy as np
from adaptive_moment_batch import setup
from compatible_moment_rhs import rhs

class CompatibleMomentTests(unittest.TestCase):
    def test_rest_pressure_has_zero_covariance_tangent(self):
        for a in (.001,.1):
            m,y=setup('wave',3,a);dy,d=rhs(0.,y,m)
            Pi=dy[:-4].reshape(9,20)[:,15:19]
            self.assertEqual(float(np.max(abs(Pi))),0.)
            self.assertLess(d['mass_flux'],3e-13)
    def test_donor_mixing_covariance_tangent_nonnegative(self):
        m,y=setup('wave',3,.01);a=y[:-4].reshape(9,20);raw=np.column_stack([m,a[:,3:]])
        vel=np.random.default_rng(732).normal(0,.03,(9,2))
        for i,u in enumerate(vel):
            raw[i,7:9]=m[i]*u;raw[i,9:13]=np.outer(raw[i,1:3],u).ravel()
            raw[i,13:17]=(m[i]*np.outer(u,u)).ravel();raw[i,17]+=.5*m[i]*(u@u)
        a[:,3:]=raw[:,1:]
        dy,d=rhs(0.,y,m);dr=dy[:-4].reshape(9,20)[:,3:]
        for i,u in enumerate(vel):
            pdot=dr[i,6:8];pidot=dr[i,12:16].reshape(2,2)
            Rdot=pidot-np.outer(pdot,u)-np.outer(u,pdot)
            self.assertGreaterEqual(np.linalg.eigvalsh((Rdot+Rdot.T)/2).min()/m[i],-1e-12)
    def test_global_pressure_work_ledger(self):
        m,y=setup('wave',3,.1);dy,d=rhs(0.,y,m);dr=dy[:-4].reshape(9,20)[:,3:]
        totals=np.r_[dr[:,6:8].sum(0),dr[:,16].sum(),(dr[:,9]-dr[:,10]).sum()]
        np.testing.assert_allclose(totals,dy[-4:],atol=1e-14,rtol=0)
if __name__=='__main__':unittest.main()

import unittest
import numpy as np
from adaptive_moment_batch import setup
from compatible_moment_rhs import rhs
from factor_moment_evolve import rebuild,advance_factors

class FactorMomentTests(unittest.TestCase):
    def test_covariance_equation_matches_raw_moments(self):
        for kind in ('wave','shear'):
            mass,y=setup(kind,3,.01);L=[np.sqrt(m)*np.array([[.01,0.],[.003,.007]]) for m in mass]
            y=rebuild(y,mass,L)
            dy,d=rhs(0.,y,mass,linear_pressure=True,covariance_factors=L,return_covariance_terms=True)
            a=y[:-4].reshape(9,20);da=dy[:-4].reshape(9,20)
            for i,m in enumerate(mass):
                q=a[i,3:5];Q=a[i,5:9].reshape(2,2);p=a[i,9:11];D=a[i,11:15].reshape(2,2)
                qd=da[i,3:5];Qd=da[i,5:9].reshape(2,2);pd=da[i,9:11];Dd=da[i,11:15].reshape(2,2);Pid=da[i,15:19].reshape(2,2)
                S=Q-np.outer(q,q)/m;C=D-np.outer(q,p)/m;B=np.linalg.solve(S,C).T
                Sd=Qd-(np.outer(qd,q)+np.outer(q,qd))/m;Cd=Dd-(np.outer(qd,p)+np.outer(q,pd))/m
                Bd=np.linalg.solve(S,(Cd.T-B@Sd).T).T
                Td=Pid-(np.outer(pd,p)+np.outer(p,pd))/m
                Rd=Td-Bd@S@B.T-B@Sd@B.T-B@S@Bd.T
                A=d['cov_A'][i];F=d['cov_Qfactor'][i];R=L[i]@L[i].T
                np.testing.assert_allclose(Rd,A@R+R@A.T+F@F.T,atol=1e-13,rtol=0)
    def test_factor_update_is_gram(self):
        L=[np.zeros((2,2))];d={'cov_A':[np.array([[.2,3.],[-2.,-.1]])],'cov_Qfactor':[np.array([[1.,.2],[.3,1.]])]}
        out=advance_factors(L,d,.1)[0]
        self.assertGreater(np.linalg.eigvalsh(out@out.T).min(),0.)
if __name__=='__main__':unittest.main()

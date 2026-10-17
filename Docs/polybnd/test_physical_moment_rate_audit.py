import unittest
import numpy as np
from physical_moment_rate_audit import centered_rates


class CenteredRateAuditTests(unittest.TestCase):
    def test_stationary_reference_adapter_supplies_time_derivatives(self):
        from vortex_rhs_consistency import time_fields
        values=time_fields(np.array([[.42,.51],[.2,.7]]),0.,.5,True)
        self.assertEqual(len(values),6)
        for value in values[3:]:np.testing.assert_array_equal(value,np.zeros_like(value))

    def fixture(self):
        m=1.3;c=np.array([.4,.6]);v=np.array([.2,-.1]);S=np.array([[.02,.003],[.003,.04]]);C=np.array([[.003,.002],[-.004,.006]])
        Pi=m*np.outer(v,v)+C.T@np.linalg.solve(S,C)+np.diag([.001,.003])
        return np.r_[m,m*c,(S+m*np.outer(c,c)).ravel(),m*v,(C+m*np.outer(c,v)).ravel(),Pi.ravel(),2.]

    def test_centered_derivatives_match_finite_differences_including_variable_mass(self):
        r=self.fixture();d=np.random.default_rng(718).normal(size=18)*.01;d[0]=.02
        def values(x):
            m=x[0];c=x[1:3]/m;v=x[7:9]/m;S=x[3:7].reshape(2,2)-m*np.outer(c,c);C=x[9:13].reshape(2,2)-m*np.outer(c,v)
            return dict(centroid=c,velocity=v,density_covariance=S/m,velocity_gradient=np.linalg.solve(S,C).T,residual_velocity_covariance=(x[13:17].reshape(2,2)-m*np.outer(v,v)-C.T@np.linalg.solve(S,C))/m)
        h=1e-5;lo=values(r-h*d);hi=values(r+h*d)
        for k,v in centered_rates(r,d).items():np.testing.assert_allclose(v,(hi[k]-lo[k])/(2*h),atol=3e-10,rtol=3e-8)

    def test_uniform_rescaling_does_not_change_specific_moments(self):
        r=self.fixture()
        for v in centered_rates(r,.3*r).values():np.testing.assert_allclose(v,0.,atol=2e-15)


if __name__=='__main__':unittest.main()

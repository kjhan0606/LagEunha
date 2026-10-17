import unittest
from unittest.mock import patch
import numpy as np
from adaptive_moment_batch import setup
from periodic2d_hllc import quadrature
from factor_moment_evolve import rebuild
from quadratic_velocity_reconstruction import prepare,velocity_samples,covariance_budget,projected_budget
from kinetic_factor_rhs import ALPHA,rhs
from smooth_mass_path import evaluate_smooth_path
import moment_ale_evolve as model


def exact_velocity(points):
    x,y=np.asarray(points).T
    return np.column_stack([.2+.1*x+.3*x*y+.2*y*y,-.1+.2*y-.15*x*x+.1*x*y])


def fixture(covariance_scale=1.,field=exact_velocity):
    mass,state=setup('wave',3,0.);a=state[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],0.);factors=[]
    for i,poly in enumerate(cells):
        points,w=quadrature(poly);v=field(points);m=w.sum();q=w@points;Q=points.T@(w[:,None]*points)
        p=w@v;D=points.T@(w[:,None]*v);Pi=v.T@(w[:,None]*v)
        S=Q-np.outer(q,q)/m;C=D-np.outer(q,p)/m;B=np.linalg.solve(S,C).T
        residual=v-p/m-(points-q/m)@B.T;_,R=np.linalg.qr(np.sqrt(w)[:,None]*residual,mode='reduced')
        factors.append(covariance_scale*R.T);mass[i]=m
        raw=np.r_[m,q,Q.ravel(),p,D.ravel(),Pi.ravel(),.6*m/(model.GAMMA-1)+.5*np.trace(Pi)]
        a[i,3:]=raw[1:]
    state=rebuild(state,mass,factors);a=state[:-4].reshape(-1,20)
    states=[model.decode(poly,np.r_[m,row[3:]]) for poly,m,row in zip(cells,mass,a)]
    return mass,state,factors,cells,edges,states


class QuadraticVelocityTests(unittest.TestCase):
    def test_general_quadratic_mean_improves_and_all_kinetic_moments(self):
        mass,state,L,cells,edges,states=fixture();before=state.copy()
        diag=prepare(states,cells,edges,mass,L)
        self.assertEqual(diag['velocity_quadratic_cells'],9)
        self.assertLess(diag['velocity_quadratic_covariance_split_error'],1e-10)
        for s,f,m in zip(states,L,mass):
            points=s['density'].points;dm=s['density'].weighted_mass
            reconstructed=np.asarray([model.velocity(s,p) for p in points])
            affine=s['v']+(points-s['c'])@s['B'].T
            error=lambda v:float(np.sum(dm[:,None]*(v-exact_velocity(points))**2))
            self.assertLess(error(reconstructed),error(affine))
            samples,_=velocity_samples(s,f,m,points,ALPHA);mean=samples.mean(1)
            p=dm@mean;D=points.T@(dm[:,None]*mean)
            Pi=np.einsum('n,nki,nkj->ij',dm/4,samples,samples)
            np.testing.assert_allclose(p,s['raw'][7:9],atol=2e-12,rtol=0.)
            np.testing.assert_allclose(D,s['raw'][9:13].reshape(2,2),atol=2e-12,rtol=0.)
            np.testing.assert_allclose(Pi,s['raw'][13:17].reshape(2,2),atol=2e-12,rtol=0.)
        np.testing.assert_array_equal(state,before)

    def test_regularized_projected_mean_is_linear_in_new_tiny_covariance(self):
        candidate=np.array([[1.,0.,0.],[0.,0.,0.]])
        ratios=[]
        for epsilon in [1e-8,1e-10,1e-12]:
            L=np.diag([np.sqrt(epsilon),0.]);H,W,T=projected_budget(L,1.,candidate,np.eye(3),regularized=True)
            self.assertLess(np.linalg.norm(H),65**2*epsilon)
            ratios.append(np.linalg.norm(H)/epsilon)
            np.testing.assert_allclose(H@H.T+L@T@T.T@L.T,L@L.T,atol=1e-20,rtol=1e-12)
        self.assertLess(abs(ratios[-1]/ratios[-2]-1),.001)

    def test_regularized_projection_still_resolves_a_full_quadratic_budget(self):
        H0=np.array([[.6,.1,.2],[.1,.5,-.1]]);L=np.linalg.cholesky(H0@H0.T)
        H,W,T=projected_budget(L,1.,H0,np.eye(3),regularized=True)
        np.testing.assert_allclose(H,H0,atol=2e-11,rtol=0.)
        np.testing.assert_allclose(H@H.T+L@T@T.T@L.T,L@L.T,atol=2e-14,rtol=0.)

    def test_regularized_projection_keeps_the_full_covariance_identity(self):
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_BUDGET':'projected_regularized'}):
            self.test_mass_only_audit_matches_full_quadratic_rhs()
            self.test_positive_factor_equation_matches_raw_covariance_derivative()

    def test_projected_budget_resolves_both_quadratic_directions(self):
        mass,state,L,cells,edges,states=fixture()
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_BUDGET':'projected'}):prepare(states,cells,edges,mass,L)
        for s,f,m in zip(states,L,mass):
            points=s['density'].points;dm=s['density'].weighted_mass
            reconstructed=np.asarray([model.velocity(s,p) for p in points])
            np.testing.assert_allclose(reconstructed,exact_velocity(points),atol=2e-8,rtol=0.)
            samples,_=velocity_samples(s,f,m,points,ALPHA)
            np.testing.assert_allclose(np.einsum('n,nki,nkj->ij',dm/4,samples,samples),s['raw'][13:17].reshape(2,2),atol=2e-12,rtol=0.)

    def test_projected_rank_limit_and_factor_gauge(self):
        candidate=np.array([[.6,.2,.1],[0.,0.,0.]])
        H0,_,_=projected_budget(np.diag([1.,0.]),1.,candidate,np.eye(3))
        rotation=np.array([[.8,-.6],[.6,.8]])
        for epsilon in [1e-3,1e-5,1e-7]:
            L=np.diag([1.,epsilon**1.5]);current=candidate.copy();current[1]=epsilon*np.array([.6,.3,-.2])
            H,W,T=projected_budget(L,1.,current,np.eye(3));h,w,t=projected_budget(L@rotation,1.,current,np.eye(3))
            self.assertLess(np.linalg.norm(H-H0),100*epsilon**3+1e-14)
            np.testing.assert_allclose(H@H.T+L@T@T.T@L.T,L@L.T,atol=3e-15,rtol=0.)
            np.testing.assert_allclose(h,H,atol=3e-15,rtol=0.)
            np.testing.assert_allclose(L@rotation@t,L@T,atol=3e-15,rtol=0.)

    def test_projected_closure_full_rhs_matches_independent_geometry(self):
        mass,state,L,cells,edges,states=fixture()
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_BUDGET':'projected'}):
            self.test_mass_only_audit_matches_full_quadratic_rhs()
            self.test_positive_factor_equation_matches_raw_covariance_derivative()

    def test_rank_one_quadratic_mean_is_reproduced(self):
        def field(points):
            x,y=np.asarray(points).T
            v=.2+.1*x+.3*x*y+.2*y*y
            return np.column_stack([v,.3*v-.1])
        mass,state,L,cells,edges,states=fixture(field=field)
        prepare(states,cells,edges,mass,L)
        for s in states:
            points=s['density'].points
            reconstructed=np.asarray([model.velocity(s,p) for p in points])
            np.testing.assert_allclose(reconstructed,field(points),atol=2e-9,rtol=0.)

    def test_rank_transition_preserves_longitudinal_curvature(self):
        candidate=np.array([[.6,.2,.1],[0.,0.,0.]])
        H0,_,_=covariance_budget(np.diag([1.,0.]),1.,candidate,np.eye(3))
        for epsilon in [1e-3,1e-5,1e-7]:
            L=np.diag([1.,epsilon**1.5]);current=candidate.copy()
            current[1]=epsilon*np.array([.6,.3,-.2])
            H,W,T=covariance_budget(L,1.,current,np.eye(3))
            self.assertLess(np.linalg.norm(H-H0),2*epsilon**3+2e-15)
            np.testing.assert_allclose(H@H.T+L@T@T.T@L.T,L@L.T,atol=2e-15,rtol=0.)
            np.testing.assert_allclose(L@T,(L@T).T,atol=2e-15,rtol=0.)

    def test_isotropic_limit_and_factor_gauge_are_continuous(self):
        candidate=np.array([[.6,.2,.1],[.4,-.1,.2]])
        angle=.7;U=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        for epsilon in [1e-3,1e-5,1e-7,0.]:
            L=U@np.diag([1.,1-epsilon])
            H,W,T=covariance_budget(L,1.,candidate,np.eye(3))
            self.assertLess(np.linalg.norm(H),3*epsilon+1e-15)
            other=L@U
            h,w,t=covariance_budget(other,1.,candidate,np.eye(3))
            np.testing.assert_allclose(h,H,atol=2e-15,rtol=0.)
            np.testing.assert_allclose(other@t,L@T,atol=2e-15,rtol=0.)

    def test_curvature_is_constrained_by_unchanged_covariance(self):
        mass,state,L,cells,edges,states=fixture(.03);before=state.copy()
        diag=prepare(states,cells,edges,mass,L)
        self.assertLessEqual(diag['velocity_quadratic_max_resolved_fraction'],1+1e-10)
        self.assertLess(diag['velocity_quadratic_covariance_split_error'],1e-10)
        for s,f in zip(states,L):
            residual=f@s['velocity_quadratic']['residual_map']
            self.assertGreaterEqual(np.linalg.eigvalsh(residual@residual.T).min(),-1e-20)
        np.testing.assert_array_equal(state,before)

    def test_nearly_rank_one_primary_factors_keep_the_budget(self):
        mass,state,L,cells,edges,states=fixture()
        L=[np.sqrt(m)*np.array([[.01,0.],[.003,1e-14]]) for m in mass]
        state=rebuild(state,mass,L);a=state[:-4].reshape(-1,20)
        states=[model.decode(p,np.r_[m,row[3:]]) for p,m,row in zip(cells,mass,a)]
        before=state.copy();diag=prepare(states,cells,edges,mass,L)
        self.assertLess(diag['velocity_quadratic_covariance_split_error'],1e-9)
        np.testing.assert_array_equal(state,before)

    def test_actual_failed_tiny_covariance_fixture(self):
        L=np.array([[-3.7043669600309045e-14,0.],[-7.437354445812192e-15,-1.0418638394587233e-23]])
        candidate=np.array([[2.3236700035987596e-8,7.037267454129931e-9,-3.4322279647129313e-9],
            [-4.689483694407818e-10,-1.4817299986258976e-10,6.78573023430591e-11]])
        G=np.array([[.0711995750523997,-.0020860940589247034,-7.267870327929381e-5],
            [0.,-.08324514982034685,-.0020519852003535576],[0.,0.,.07814755160270612]])
        mass=.0032744613402471206;H,W,T=covariance_budget(L,mass,candidate,G)
        R=L@L.T;reconstructed=mass*H@(G.T@G)@H.T+L@T@T.T@L.T
        self.assertLess(np.linalg.norm(reconstructed-R)/np.linalg.norm(R),1e-12)
        self.assertLess(np.max(np.linalg.svd(W@G.T,compute_uv=False)),1.)

    def test_mass_only_audit_matches_full_quadratic_rhs(self):
        mass,state,L,cells,edges,states=fixture()
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_VELOCITY':'1','LAGFORCE_VECTOR_KINETIC_ROOT':'1'}):
            dy,diag=rhs(0.,state,mass,covariance_factors=L)
            rates=dy[:-4].reshape(-1,20);flux,audit=evaluate_smooth_path(0.,state,L,mass,(rates[:,:2],rates[:,2]))
        np.testing.assert_allclose(flux,diag['mass_rates'],atol=2e-12,rtol=0.)
        self.assertLess(audit['generator_defect'],1e-12)

    def test_positive_factor_equation_matches_raw_covariance_derivative(self):
        mass,state,L,cells,edges,states=fixture()
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_VELOCITY':'1','LAGFORCE_VECTOR_KINETIC_ROOT':'1'}):
            dy,diag=rhs(0.,state,mass,covariance_factors=L)
        a=state[:-4].reshape(-1,20);da=dy[:-4].reshape(-1,20)
        for i,m in enumerate(mass):
            q=a[i,3:5];Q=a[i,5:9].reshape(2,2);p=a[i,9:11];D=a[i,11:15].reshape(2,2)
            qd=da[i,3:5];Qd=da[i,5:9].reshape(2,2);pd=da[i,9:11];Dd=da[i,11:15].reshape(2,2);Pid=da[i,15:19].reshape(2,2)
            S=Q-np.outer(q,q)/m;C=D-np.outer(q,p)/m;B=np.linalg.solve(S,C).T
            Sd=Qd-(np.outer(qd,q)+np.outer(q,qd))/m;Cd=Dd-(np.outer(qd,p)+np.outer(q,pd))/m
            Bd=np.linalg.solve(S,(Cd.T-B@Sd).T).T;Td=Pid-(np.outer(pd,p)+np.outer(p,pd))/m
            Rd=Td-Bd@S@B.T-B@Sd@B.T-B@S@Bd.T
            R=L[i]@L[i].T;A=diag['cov_A'][i];F=diag['cov_Qfactor'][i]
            expected=A@R+R@A.T-L[i]@diag['cov_K'][i]@L[i].T+F@F.T
            np.testing.assert_allclose(Rd,expected,atol=5e-11,rtol=0.)
            self.assertGreaterEqual(np.linalg.eigvalsh(diag['cov_K'][i]).min(),-1e-12)

    def test_zero_covariance_affine_patch_is_preserved(self):
        mass,state=setup('shear',3,0.);L=[np.zeros((2,2)) for _ in mass]
        state=rebuild(state,mass,L)
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_VELOCITY':'0','LAGFORCE_VECTOR_KINETIC_ROOT':'1'}):old,_=rhs(0.,state,mass,covariance_factors=L)
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_VELOCITY':'1','LAGFORCE_VECTOR_KINETIC_ROOT':'1'}):new,diag=rhs(0.,state,mass,covariance_factors=L)
        np.testing.assert_allclose(new,old,atol=1e-12,rtol=0.)
        self.assertEqual(diag['velocity_quadratic_cells'],0)

    def test_unsupported_linear_bin_transport_is_rejected(self):
        mass,state,L,*_=fixture()
        with patch.dict('os.environ',{'LAGFORCE_QUADRATIC_VELOCITY':'1'}):
            with self.assertRaisesRegex(ValueError,'continuous moment density'):
                rhs(0.,state,mass,covariance_factors=L,exact_bin_transport=True)


if __name__=='__main__':unittest.main()

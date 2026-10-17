"""Run with explicit actual-C geometry library on a compute host."""
import argparse,json
from pathlib import Path
import numpy as np
from power_cells_3d import Geometry3D
from moment_density_3d import quadrature
from moment_gravity_kick import moments
from moment_state_3d import pack,unpack,ALPHA
from kinetic_moment_rhs_3d import rhs


def initial(geometry,sites,weights,velocity,dispersion=0.,bounds=None):
    polys,_,_=geometry.build(sites,weights,bounds);raw=[];factors=[]
    for poly in polys:
        x,w=quadrature(poly);v=velocity(x)[:,None,:]+dispersion*ALPHA[None,:,:]
        r=moments(x,w,v,1.5*w.sum());raw.append(pack(*(r[k] for k in ('m','q','Q','p','D','Pi','E'))));factors.append(r['residual_factor'])
    return np.array(raw),np.array(factors)


def covariance_error(raw,dr,L,diag):
    error=0.
    for i,r in enumerate(raw):
        m,q,Q,p,D,Pi,E=unpack(r);_,qd,Qd,pd,Dd,Pid,_=unpack(dr[i])
        S=Q-np.outer(q,q)/m;C=D-np.outer(q,p)/m;B=np.linalg.solve(S,C).T
        Sd=Qd-(np.outer(qd,q)+np.outer(q,qd))/m;Cd=Dd-(np.outer(qd,p)+np.outer(q,pd))/m
        Bd=np.linalg.solve(S,(Cd.T-B@Sd).T).T;Td=Pid-(np.outer(pd,p)+np.outer(p,pd))/m
        Rd=Td-Bd@S@B.T-B@Sd@B.T-B@S@Bd.T
        A=diag['cov_A'][i];F=diag['cov_Qfactor'][i];R=L[i]@L[i].T
        predicted=A@R+R@A.T-L[i]@diag['cov_K'][i]@L[i].T+F@F.T
        error=max(error,float(np.max(abs(Rd-predicted))))
    return error


def run(library):
    geometry=Geometry3D(library);rng=np.random.default_rng(271)
    sites=np.array([[i,j,k] for i in (.25,.75) for j in (.25,.75) for k in (.25,.75)])+rng.uniform(-.035,.035,(8,3));weights=np.zeros(8);rows=[]
    for velocity in (np.zeros(3),np.array([.2,-.1,.3])):
        for t in (0.,.1):
            x=sites+t*velocity;bounds=np.array([[0.,0.,0.],[1.,1.,1.]])+t*velocity
            raw,L=initial(geometry,x,weights,lambda p:np.tile(velocity,(len(p),1)),bounds=bounds)
            dr,gen,z,diag=rhs(geometry,x,weights,raw,L,bounds,lambda p:velocity,face_order=3)
            exact=np.zeros_like(raw)
            for i,r in enumerate(raw):
                m,q,Q,p,D,Pi,E=unpack(r);exact[i,1:4]=p;exact[i,4:13]=(np.outer(q,velocity)+np.outer(velocity,q)).ravel();exact[i,16:25]=Pi.ravel()
            error=float(np.max(abs(dr-exact)));assert error<2e-11
            np.testing.assert_allclose(gen,np.tile(velocity,(len(x),1)),atol=2e-13,rtol=0);assert max(abs(z))<2e-12
            rows.append(dict(case='uniform_translation',t=t,velocity=velocity.tolist(),max_raw_rhs_error=error,mass_flux=diag['mass_flux'],partition_error=diag['partition_error']))
    def nonuniform(p):return .02*np.column_stack([np.sin(np.pi*p[:,0])*np.cos(np.pi*p[:,1]),np.sin(np.pi*p[:,1])*np.cos(np.pi*p[:,2]),np.sin(np.pi*p[:,2])*np.cos(np.pi*p[:,0])])
    raw,L=initial(geometry,sites,weights,nonuniform,dispersion=.003)
    dr,gen,z,diag=rhs(geometry,sites,weights,raw,L,face_order=3)
    totals=np.r_[dr[:,13:16].sum(axis=0),dr[:,34].sum(),(dr[:,[21,22,17]]-dr[:,[23,18,19]]).sum(axis=0)]
    # D indices: yz=21, zy=23; zx=22, xz=18; xy=17, yx=19.
    ledger_error=float(max(abs(totals-diag['ledger'])));cov_error=covariance_error(raw,dr,L,diag)
    assert ledger_error<2e-12 and cov_error<2e-12 and diag['mass_flux']<2e-13
    rows.append(dict(case='nonuniform',ledger_error=ledger_error,covariance_identity_error=cov_error,mass_flux=diag['mass_flux'],weight_rate_max=float(max(abs(z))),generator_velocity_defect=diag['generator_velocity_defect']))
    return dict(passed=True,rows=rows,scope='Actual 3D C geometry plus density/moment/kinetic RHS; cell-constant pressure, no time integration yet')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('library');args=parser.parse_args();result=run(args.library)
    Path('coupled_rhs_3d_report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)

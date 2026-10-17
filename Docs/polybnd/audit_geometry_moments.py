"""Static GCL/moment audits of saved shear states and an initial RHS tangent.
No evolution or correction. The affine reference is only a diagnostic.
"""
import json
from pathlib import Path
import numpy as np
from adaptive_moment_batch import setup
from periodic2d_hllc import quadrature
import moment_ale_evolve as model
from compatible_moment_rhs import rhs as donor_rhs
from kinetic_factor_rhs import rhs as kinetic_rhs

GROUPS={'mass':slice(0,1),'q':slice(1,3),'Q':slice(3,7),'p':slice(7,9),'D':slice(9,13),'Pi':slice(13,17),'E':slice(17,18)}

def exact_raw(polys):
    rows=[]
    for poly in polys:
        r,w=quadrature(poly);v=np.column_stack([.2+1.1*r[:,1],np.zeros(len(w))]);Pi=v.T@(w[:,None]*v)
        rows.append(np.r_[sum(w),w@r,(r.T@(w[:,None]*r)).ravel(),w@v,(r.T@(w[:,None]*v)).ravel(),Pi.ravel(),sum(w)*.6/(model.GAMMA-1)+.5*np.trace(Pi)])
    return np.array(rows)

def compare_saved(path):
    mass,_=setup('shear',4,0.);z=json.loads(path.read_text());t=z['t'];mass=np.array(z['mass']);a=np.array(z['state'])[:-4].reshape(-1,20)
    polys,_=model.geometry(a[:,:2],a[:,2],t);exact=exact_raw(polys);raw=np.column_stack([mass,a[:,3:]])
    reconstructed=[model.decode(poly,row) for poly,row in zip(polys,raw)]
    rho_error=max(float(np.max(abs(s['density'].value(s['density'].points)-1))) for s in reconstructed)
    return dict(path=str(path),time=t,max_relative_volume_mass_mismatch=float(max(abs(exact[:,0]-mass)/mass)),
                moment_max_abs={k:float(np.max(abs(raw[:,sl]-exact[:,sl]))) for k,sl in GROUPS.items()},
                density_departure_from_one=rho_error,pressure_departure=max(abs(s['P']-.6) for s in reconstructed),
                fitted_moment_residual=max(s['density'].residual for s in reconstructed))

def tangent():
    mass,y=setup('shear',3,0.);N=len(mass);L=[np.zeros((2,2)) for _ in mass];answer={}
    for name,rhs in [('donor',donor_rhs),('kinetic',kinetic_rhs)]:
        dy,diag=rhs(0.,y,mass,linear_pressure=True,covariance_factors=L,return_covariance_terms=True)
        a=y[:-4].reshape(N,20);da=dy[:-4].reshape(N,20);rates=np.column_stack([np.zeros(N),da[:,3:]])
        errors=[]
        for h in (1e-4,5e-5,2.5e-5):
            minus,_=model.geometry(a[:,:2]-h*da[:,:2],a[:,2]-h*da[:,2],-h)
            plus,_=model.geometry(a[:,:2]+h*da[:,:2],a[:,2]+h*da[:,2],h)
            fd=(exact_raw(plus)-exact_raw(minus))/(2*h)
            errors.append(dict(h=h,groups={k:float(np.max(abs(fd[:,sl]-rates[:,sl]))) for k,sl in GROUPS.items()}))
        answer[name]=dict(mass_flux=diag['mass_flux'],centered_geometry_derivative_errors=errors)
    return answer

def main():
    folder=Path('/gpfs/kjhan/LagForce/Codex_B_factor_moments_20261010_002319')
    report=dict(saved=[compare_saved(folder/f'shear_a0_dt{dt:g}.checkpoint.json') for dt in (.001,.0005)],tangent=tangent())
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()

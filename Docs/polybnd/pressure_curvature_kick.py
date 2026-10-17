"""Local Euler acceleration test at rest under a positive quadratic pressure.

The interior exact initial acceleration is affine, with gradient -H. No exact
acceleration is supplied to the RHS. Outer-wall compatibility is deliberately
not used to infer a global smooth evolution from this instantaneous test.
"""
import json,time,os
from pathlib import Path
import numpy as np
from adaptive_moment_batch import setup
from adaptive_coupled_mass import atomic
from factor_moment_evolve import rebuild
from periodic2d_hllc import quadrature
from simple_wave_evolution import install_domain
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic

HESSIAN=np.array([[.1,.03],[.03,-.04]])


def fixture(n=8,regular=False):
    install_domain();mass,y=setup('wave',n,0.);a=y[:-4].reshape(-1,20)
    if regular:
        xx,yy=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n);a[:,:2]=np.column_stack([xx.ravel(),yy.ravel()])
    cells,_=model.geometry(a[:,:2],a[:,2],0.)
    for i,poly in enumerate(cells):
        points,dm=quadrature(poly);z=points-.5;P=.8+.5*np.einsum('ni,ij,nj->n',z,HESSIAN,z)
        mass[i]=dm.sum();a[i,3:5]=dm@points;a[i,5:9]=(points.T@(dm[:,None]*points)).ravel()
        a[i,9:19]=0.;a[i,19]=dm@P/(model.GAMMA-1)
    factors=[np.zeros((2,2)) for _ in mass]
    return mass,rebuild(y,mass,factors),factors


def inspect(n=8,regular=False):
    mass,y,L=fixture(n,regular);dy,diag=kinetic.rhs(0.,y,mass,covariance_factors=L)
    a=y[:-4].reshape(-1,20);da=dy[:-4].reshape(-1,20);rows=[]
    for i,m in enumerate(mass):
        c=a[i,3:5]/m
        if np.any(c<.2) or np.any(c>.8):continue
        S=a[i,5:9].reshape(2,2)-m*np.outer(c,c)
        Cdot=da[i,11:15].reshape(2,2)-np.outer(c,da[i,9:11])
        gradient=np.linalg.solve(S,Cdot).T;acceleration=da[i,9:11]/m
        rows.append(dict(cell=i,center=c,mean_acceleration_error=float(np.max(abs(acceleration+HESSIAN@(c-.5)))),
            affine_acceleration_gradient=gradient,gradient_error=float(np.max(abs(gradient+HESSIAN)))))
    return dict(n=n,regular_sites=regular,pressure_reconstruction=diag.get('pressure_reconstruction','linear'),cells_tested=len(rows),
        maximum_mean_acceleration_error=max(v['mean_acceleration_error'] for v in rows),maximum_gradient_error=max(v['gradient_error'] for v in rows),
        exact_affine_acceleration_gradient=-HESSIAN,rows=rows,mass_residual=diag['mass_flux'],
        generator_defect=diag['generator_velocity_defect'],
        note='Instantaneous interior Euler source test, rho=1,u=0,P=.8+.5*(x-.5)^T H*(x-.5). Pressure initialized through cell total energy only. Exact acceleration enters scoring only. No global wall-compatible trajectory or physical accuracy pass is asserted.')


def rotating_patch(n=8,omega=.2):
    """Full moment Reynolds derivative for a locally stationary rotating flow."""
    from simple_wave_exact import polygon_rule
    from simple_wave_rhs_consistency import exact_rates
    mass,y,L=fixture(n);a=y[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],0.)
    def fields(points,t,amplitude,with_time_derivatives=False):
        z=np.asarray(points)-.5;v=amplitude*np.column_stack([-z[:,1],z[:,0]]);P=.8+.5*amplitude**2*np.sum(z*z,axis=1);rho=np.ones(len(z))
        if with_time_derivatives:return rho,v,P,np.zeros_like(rho),np.zeros_like(v),np.zeros_like(P)
        return rho,v,P
    def raw(poly,t,amplitude,order):
        points,w=polygon_rule(poly,order);rho,v,P=fields(points,t,amplitude);Pi=v.T@(w[:,None]*v)
        r=np.r_[w.sum(),w@points,(points.T@(w[:,None]*points)).ravel(),w@v,(points.T@(w[:,None]*v)).ravel(),Pi.ravel(),w@P/(model.GAMMA-1)+.5*np.trace(Pi)]
        return r,points,w,rho,v,P
    for i,p in enumerate(cells):
        r,*_=raw(p,0.,omega,8);mass[i]=r[0];a[i,3:]=r[1:]
    y=rebuild(y,mass,L);dy,diag=kinetic.rhs(0.,y,mass,covariance_factors=L);a=y[:-4].reshape(-1,20);da=dy[:-4].reshape(-1,20)
    reference,_,_=exact_rates(cells,edges,a[:,:2],da[:,:2],da[:,2],0.,omega,8,fields,raw)
    numerical=np.column_stack([np.zeros(len(mass)),da[:,3:]])
    centers=a[:,3:5]/mass[:,None];interior=np.all((centers>.2)&(centers<.8),axis=1)
    difference=np.max(abs(numerical[interior]-reference[interior])/mass[interior,None],axis=0)
    return dict(n=n,omega=omega,pressure_reconstruction=diag.get('pressure_reconstruction','linear'),interior_cells=int(interior.sum()),
        maximum_raw_moment_rate_error_per_mass=float(difference.max()),raw_moment_rate_errors_per_mass=difference,mass_residual=diag['mass_flux'],generator_defect=diag['generator_velocity_defect'],
        note='Global affine rotating velocity and matching quadratic pressure: stationary Euler fields, moving-cell Reynolds rates. Interior only; reflecting-box wall compatibility not assumed. Exact fields enter initialization/reference only. All18 raw moment derivatives checked.')


if __name__=='__main__':
    start=time.monotonic();r=inspect();r['seconds']=time.monotonic()-start
    atomic('pressure_curvature_kick_report.json',r);print(r,flush=True)

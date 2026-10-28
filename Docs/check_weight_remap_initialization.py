"""Independent quadrature and density-change bound for the positive lift."""
import importlib.util
import json
from pathlib import Path
import numpy as np


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod


def main():
    root=Path(__file__).parent
    base=Path('/gpfs/kjhan/LagForce')
    mod=load('moments',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
    geo=load('geometry',base/'Codex_weighted_compression_geometry_20261017/audit.py')
    fit=base/'Codex_weight_remap_initialization_20261018'
    state=np.load(fit/'candidate.npz');ids=state['ids'];theta=state['theta'];scale=state['scale']
    a=np.loadtxt(base/'Codex_gradient_operator_20261017/operator/rhs.dat');a=a[np.argsort(a[:,0])]
    x=a[:,1:3]-3;vel=a[:,3:5];mass=a[ids,7]
    polys=mod.polygons(x,np.zeros(len(x)),ids,geo.clip)
    volume=np.array([mod.area(p) for p in polys]);phi0=mod.potential(x[ids]);rho=mass/volume
    matrices=[];rhs=[];means=[];norms=[];quadratures=[]
    for k,p in enumerate(polys):
        xy,dv=mod.quadrature(p,4);r=xy-x[ids[k]];v=vel[ids[k]]
        f=np.column_stack((r,r[:,0]*v[1]-r[:,1]*v[0],mod.potential(xy)-phi0[k]))/scale
        mean=np.sum(dv[:,None]*f,axis=0)/volume[k]
        g=rho[k]*(f-mean)
        matrices.append(g.T@(dv[:,None]*g));rhs.append(rho[k]*(dv@f));means.append(mean)
        norms.append(float(dv@np.exp(f@theta)/volume[k]));quadratures.append((f,dv))
    G=np.sum(matrices,axis=0);b=-np.sum(rhs,axis=0);lam=np.linalg.solve(G,b)
    bound=float(np.sqrt((b@lam)/volume.sum()))
    factors=[];l2total=0.;achieved=np.zeros(4)
    for k,(f,dv) in enumerate(quadratures):
        delta=(rho[k]*(f-means[k]))@lam;factors.extend(1+delta)
        l2total+=dv@(delta*delta);achieved+=rho[k]*(dv@(delta[:,None]*f))
    assert np.max(abs(achieved-b))<1e-11
    assert abs(np.sqrt(l2total/volume.sum())-bound)<1e-12
    # The unconstrained L2 minimizer happened to be positive at the rule's
    # points. Also bound its continuous density throughout every polygon:
    # the linear term attains extrema at vertices and spline phi is monotone
    # with radius. Separate extrema give a conservative bound on their sum.
    lower=[];upper=[];coeff=lam/scale
    for k,p in enumerate(polys):
        r=p-x[ids[k]];v=vel[ids[k]]
        linear=r@np.array([coeff[0]+coeff[2]*v[1],coeff[1]-coeff[2]*v[0]])
        closest=[];inside=True
        for aa,bb in zip(p,np.roll(p,-1,axis=0)):
            edge=bb-aa
            inside &= (edge[0]*(-aa[1])-edge[1]*(-aa[0])>=0)
            t=np.clip(-(aa@edge)/(edge@edge),0,1) if edge@edge else 0.
            closest.append(np.linalg.norm(aa+t*edge))
        rmin=0. if inside else min(closest);rmax=np.max(np.linalg.norm(p,axis=1))
        phirange=mod.potential(np.array([[rmin,0.],[rmax,0.]]))-phi0[k]
        offset=1-rho[k]*(means[k]@lam)
        lower.append(offset+rho[k]*(linear.min()+np.min(coeff[3]*phirange)))
        upper.append(offset+rho[k]*(linear.max()+np.max(coeff[3]*phirange)))
    checks=[];l2checks=[]
    for order in (4,8,12):
        residual=np.zeros(4);mass_error=[];l2=0.;minfactor=np.inf;maxfactor=0.
        l2res=np.zeros(4);l2mass=[];l2rms=0.
        for k,p in enumerate(polys):
            xy,dv=mod.quadrature(p,order);r=xy-x[ids[k]];v=vel[ids[k]]
            f=np.column_stack((r,r[:,0]*v[1]-r[:,1]*v[0],mod.potential(xy)-phi0[k]))/scale
            factor=np.exp(f@theta)/norms[k];dm=rho[k]*dv*factor
            mass_error.append(float(dm.sum()/mass[k]-1))
            residual+=dm@f;l2+=dv@((factor-1)**2)
            minfactor=min(minfactor,float(factor.min()));maxfactor=max(maxfactor,float(factor.max()))
            l2factor=1+rho[k]*((f-means[k])@lam)
            l2dm=rho[k]*dv*l2factor
            l2mass.append(float(l2dm.sum()/mass[k]-1))
            l2res+=l2dm@f;l2rms+=dv@((l2factor-1)**2)
        checks.append({'quadrature_order':order,'same_order4_density_normalization':True,
                       'cell_mass_max_relative_error':float(np.max(np.abs(mass_error))),
                       'global_centered_constraint_residual':(residual*scale).tolist(),
                       'density_rms_change':float(np.sqrt(l2/volume.sum())),
                       'sampled_density_min_factor':minfactor,'sampled_density_max_factor':maxfactor})
        l2checks.append({'quadrature_order':order,'cell_mass_max_relative_error':float(np.max(np.abs(l2mass))),
                         'global_centered_constraint_residual':(l2res*scale).tolist(),
                         'density_rms_change':float(np.sqrt(l2rms/volume.sum()))})
    result={'scope':'Checks the saved positive lift without refitting. Four constraints are global x/y first mass moment, angular momentum and external potential.',
            'quadrature_checks':checks,
            'unconstrained_L2_density_change_lower_bound':bound,
            'L2_optimal_min_factor':float(min(factors)),'L2_optimal_max_factor':float(max(factors)),
            'L2_continuous_density_lower_bound':float(min(lower)),
            'L2_continuous_density_upper_bound':float(max(upper)),
            'L2_independent_quadrature_checks':l2checks,
            'bound_scope':'Same order4 support and piecewise constant cell velocities, fixed per-cell mass. Bound minimizes volume-weighted squared density-factor change, allowing negative factors; hence it is a necessary lower bound for positive fits, not a theorem for all velocity/density reconstructions.',
            'interpretation':'Meeting global invariants by adjusting hidden density is not automatically a faithful restart. Compare density distortion and independent quadrature errors before adopting.'}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(root/'l2_candidate.npz',ids=ids,coefficient=lam,feature_scale=scale,
             cell_feature_mean=np.array(means),volume=volume,
             density_lower_bound=np.array(lower),density_upper_bound=np.array(upper))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()

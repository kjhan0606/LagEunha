"""Conditional point-checkpoint lifting using positive within-cell density.

Keep each old m,p,K,U, and fit global first mass moment, angular momentum and
external potential to the point convention. This changes subcell density and
must be measured, not silently accepted as a production restart conversion.
"""
import ctypes
import hashlib
import importlib.util
import json
from pathlib import Path
import time

import numpy as np


def main():
    root=Path(__file__).parent
    assert not (root/'result.json').exists()
    base=Path('/gpfs/kjhan/LagForce')
    parent=base/'Codex_weight_remap_moments_centered_20261017'
    spec=importlib.util.spec_from_file_location('moments',parent/'audit.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    geom=base/'Codex_weighted_compression_geometry_20261017/audit.py'
    spec=importlib.util.spec_from_file_location('geom',geom)
    geo=importlib.util.module_from_spec(spec);spec.loader.exec_module(geo)
    source=base/'Codex_gradient_operator_20261017/operator/rhs.dat'
    a=np.loadtxt(source);a=a[np.argsort(a[:,0])]
    x=a[:,1:3]-3;vel=a[:,3:5];ids=np.flatnonzero(np.linalg.norm(x,axis=1)<.7)
    poly=mod.polygons(x,np.zeros(len(x)),ids,geo.clip)
    started=time.monotonic();xy=[];dv=[];group=[]
    for k,p in enumerate(poly):
        q,w=mod.quadrature(p,4);xy.append(q);dv.append(w);group.extend([k]*len(w))
    xy=np.concatenate(xy);dv=np.concatenate(dv);group=np.array(group)
    n=len(ids);mass=a[ids,7];volume=np.bincount(group,weights=dv,minlength=n)
    raw=dv/volume[group];phi=mod.potential(xy);phi0=mod.potential(x[ids])
    rel=xy-x[ids[group]];v=vel[ids[group]]
    f=np.column_stack((rel,rel[:,0]*v[:,1]-rel[:,1]*v[:,0],phi-phi0[group]))
    scale=np.sqrt(np.sum((mass[group]*raw)[:,None]*f*f,axis=0)/mass.sum())
    f=f/scale;theta=np.zeros(4)

    def evaluate(t,with_hessian=True):
        z=f@t;maximum=np.full(n,-np.inf);np.maximum.at(maximum,group,z)
        q=raw*np.exp(z-maximum[group]);norm=np.bincount(group,weights=q,minlength=n)
        dm=mass[group]*q/norm[group]
        mean=np.column_stack([np.bincount(group,weights=dm*f[:,j],minlength=n)/mass for j in range(4)])
        residual=np.sum(mass[:,None]*mean,axis=0)/mass.sum()
        objective=float(np.sum(mass*(maximum+np.log(norm)))/mass.sum())
        if not with_hessian:return objective,residual,dm,None
        centered=f-mean[group]
        hessian=centered.T@((dm/mass.sum())[:,None]*centered)
        return objective,residual,dm,hessian

    history=[]
    for it in range(40):
        objective,residual,dm,hessian=evaluate(theta)
        history.append({'iteration':it,'max_scaled_constraint':float(np.max(abs(residual))),
                        'hessian_condition':float(np.linalg.cond(hessian))})
        if max(abs(residual))<2e-13:break
        step=np.linalg.solve(hessian,-residual);slope=residual@step
        alpha=1.
        for trial in range(40):
            trial_objective,trial_residual,_,_=evaluate(theta+alpha*step,False)
            if trial_objective<=objective+1e-4*alpha*slope or np.max(abs(trial_residual))<2e-13:break
            alpha*=.5
        else:raise RuntimeError('Density-fit line search failed')
        theta+=alpha*step
    else:raise RuntimeError('Density-fit Newton failed')
    assert np.min(dm)>0
    factor=dm/(mass[group]*raw)
    pressure=2.7777796321446218e-5;gamma=1.6666669845581055
    # Uniform thermal energy per volume preserves uniform pressure while the
    # density changes. No compensating transfer into heat is made.
    specific_thermal=pressure*dv/((gamma-1)*dm)
    data=np.ascontiguousarray(np.column_stack((dm,xy,v,specific_thermal,phi)))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    cpu=ctypes.CDLL(str(parent/'cpu.so'));gpu=ctypes.CDLL(str(parent/'cuda.so'))
    cpu.wm_set_threads(4)
    outputs=[]
    for fun in (cpu.wm_batch,gpu.wm_cuda_batch):
        fun.argtypes=[ctypes.c_size_t,arr,arr];fun.restype=ctypes.c_int
        out=np.empty((len(data),17));assert fun(len(data),data,out)==0;outputs.append(out)
    assert np.array_equal(outputs[0].view('u8'),outputs[1].view('u8'))
    state=np.zeros((n,17));np.add.at(state,group,outputs[0])
    thermal=state[:,15]-.5*(state[:,12]+state[:,14])
    constraints={'mass_max_relative_error':float(np.max(abs(state[:,0]/mass-1))),
                 'momentum_max_absolute_error':float(np.max(abs(state[:,6:8]-mass[:,None]*vel[ids]))),
                 'global_first_mass_moment_error':(state[:,1:3].sum(axis=0)-(mass[:,None]*x[ids]).sum(axis=0)).tolist(),
                 'global_angular_error':float(np.sum(state[:,9]-state[:,10])-np.sum(mass*(x[ids,0]*vel[ids,1]-x[ids,1]*vel[ids,0]))),
                 'global_point_potential_error':float(state[:,16].sum()-mass@phi0),
                 'pressure_max_relative_error':float(np.max(abs((gamma-1)*thermal/(volume*pressure)-1)))}
    assert constraints['mass_max_relative_error']<1e-11
    assert max(abs(constraints['global_angular_error']),abs(constraints['global_point_potential_error']))<1e-11
    assert constraints['pressure_max_relative_error']<1e-8
    result={'scope':'Positive density fit on old cells at a single canonical checkpoint; no dynamics or automatic restart conversion.',
            'constraints':constraints,'cpu_openmp4_cuda_bitwise_equal':True,
            'density_factor_min':float(factor.min()),'density_factor_max':float(factor.max()),
            'density_fractional_change_volume_weighted_rms':float(np.sqrt(dv@((factor-1)**2)/dv.sum())),
            'theta':theta.tolist(),'feature_scale':scale.tolist(),'iterations':history,
            'minimum_thermal_energy':float(thermal.min()),'seconds':time.monotonic()-started,
            'limitations':['Preserves global, not per-cell, angular/potential/position moments.',
                           'Density changes inside cells. Bounded smoothness, fluid accuracy and reconstruction at faces require separate verification.',
                           'The fit is to the selected positive quadrature. Continuous-integral accuracy is not proven by quadrature conservation.',
                           'The 17th potential moment still needs consistent gravity work in physical evolution.',
                           'Does not implement the approved Voronoi reduction or make a long-time stability claim.'],
            'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,parent/'audit.py',parent/'weight_remap_moments.h',Path(__file__))}}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(root/'candidate.npz',ids=ids,theta=theta,scale=scale,state=state,volume=volume)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()

"""Finite overlap remap of actual canonical cells; run on lageunha.

Donor fields are constant per original cell. No physical time advancement.
Positive overlap quadrature is an explicit realization, not just a PSD test.
The old point-state to extended-cell initialization is audited, never applied
silently to a production checkpoint.
"""
import ctypes
import hashlib
import importlib.util
import json
from pathlib import Path
import time

import numpy as np
from scipy.spatial import ConvexHull


def area(p):
    if len(p)<3: return 0.
    # Translation avoids cancellation falsely turning a shared edge into a
    # positive-area fragment when global coordinates are much larger.
    p=p-p[0]
    q=np.roll(p,-1,axis=0)
    return .5*np.sum(p[:,0]*q[:,1]-p[:,1]*q[:,0])


def polygons(x,w,ids,clip):
    hull=ConvexHull(np.column_stack((x,np.sum(x*x,axis=1)-w)))
    tri=hull.simplices[hull.equations[:,2]<-1e-12]
    edges=np.concatenate((tri[:,[0,1]],tri[:,[0,2]],tri[:,[1,2]]))
    edges=np.unique(np.sort(edges,axis=1),axis=0)
    nb=[[] for _ in x]
    for i,j in edges: nb[i].append(j);nb[j].append(i)
    result=[]
    for i in ids:
        p=np.array([[-3.,-3.],[3.,-3.],[3.,3.],[-3.,3.]])-x[i]
        for j in sorted(nb[i]):
            d=x[j]-x[i]
            p=clip(p,2*d,d@d+w[i]-w[j])
            assert len(p)>=3, ('empty cell',int(i))
        assert area(p)>0
        result.append(p+x[i])
    return result


def quadrature(p,order):
    # Positive tensor Gauss rule mapped onto a centroid fan of triangles.
    z,g=np.polynomial.legendre.leggauss(order);z=(z+1)/2;g=g/2
    u,v=np.meshgrid(z,z,indexing='ij');gu,gv=np.meshgrid(g,g,indexing='ij')
    u,v=u.ravel(),v.ravel();weight=(gu*gv).ravel()*(1-u)
    center=np.mean(p,axis=0);points=[];weights=[]
    for a,b in zip(p,np.roll(p,-1,axis=0)):
        da,db=a-center,b-center
        determinant=da[0]*db[1]-da[1]*db[0]
        if determinant<=0: continue
        points.append(center+u[:,None]*da+(1-u[:,None])*v[:,None]*db)
        weights.append(weight*determinant)
    return np.concatenate(points),np.concatenate(weights)


def potential(pos):
    u=np.linalg.norm(pos,axis=1)/.5
    out=np.empty(len(u));low=u<.5;mid=(u>=.5)&(u<1);high=u>=1
    z=u[low];out[low]=(16/3*z*z-48/5*z**4+32/5*z**5-14/5)/.5
    z=u[mid];out[mid]=(1/(15*z)+32/3*z*z-16*z**3+48/5*z**4-32/15*z**5-16/5)/.5
    out[high]=-1/(.5*u[high])
    return out


def main():
    root=Path(__file__).parent
    assert not (root/'result.json').exists(), 'Preserve previous result'
    start=time.monotonic();base=Path('/gpfs/kjhan/LagForce')
    source=base/'Codex_gradient_operator_20261017/operator/rhs.dat'
    code=base/'Codex_weighted_compression_geometry_20261017/audit.py'
    spec=importlib.util.spec_from_file_location('geometry',code)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    a=np.loadtxt(source);a=a[np.argsort(a[:,0])]
    x=a[:,1:3]-3;vel=a[:,3:5];ids=np.flatnonzero(np.linalg.norm(x,axis=1)<.7)
    target=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights']
    old=polygons(x,np.zeros(len(x)),ids,mod.clip)
    vol=np.array([area(p) for p in old]);mass=a[ids,7];rho=mass/vol
    assert np.max(abs(vol/a[ids,6]-1))<1e-7
    gamma=1.6666669845581055;pressure=2.7777796321446218e-5
    thermal=pressure/((gamma-1)*rho)
    cpu=ctypes.CDLL(str(root/'cpu.so'));gpu=ctypes.CDLL(str(root/'cuda.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for fun in (cpu.wm_batch,gpu.wm_cuda_batch):
        fun.argtypes=[ctypes.c_size_t,arr,arr];fun.restype=ctypes.c_int

    def call(fun,data):
        data=np.ascontiguousarray(data);out=np.empty((len(data),17))
        assert fun(len(data),data,out)==0
        assert np.isfinite(out).all()
        return out

    def atoms(p,k,order):
        xy,dv=quadrature(p,order)
        return np.column_stack((rho[k]*dv,xy,np.tile(vel[ids[k]],(len(xy),1)),
                                np.full(len(xy),thermal[k]),potential(xy)))

    def independent_old(order):
        return np.array([call(cpu.wm_batch,atoms(p,k,order)).sum(axis=0) for k,p in enumerate(old)])

    assert cpu.wm_set_threads(1)==1
    ref=independent_old(4);ref8=independent_old(8)
    lo=np.array([p.min(axis=0) for p in old]);hi=np.array([p.max(axis=0) for p in old])
    reports=[]
    for fraction in (0.,.25,1.):
        new=polygons(x,fraction*target,ids,mod.clip)
        newvol=np.array([area(p) for p in new])
        inputs=[];donors=[];receivers=[];pieces=0;coverold=np.zeros(len(ids));covernew=np.zeros(len(ids))
        phi8=0.
        for dest,p in enumerate(new):
            pmin,pmax=p.min(axis=0),p.max(axis=0)
            candidates=np.flatnonzero(np.all(hi>=pmin,axis=1)&np.all(lo<=pmax,axis=1))
            for src in candidates:
                overlap=old[src].copy()
                for aa,bb in zip(p,np.roll(p,-1,axis=0)):
                    edge=bb-aa;normal=np.array([edge[1],-edge[0]])
                    overlap=mod.clip(overlap,normal,normal@aa)
                    if len(overlap)<3: break
                av=area(overlap)
                if av<=1e-18: continue
                coverold[src]+=av;covernew[dest]+=av;pieces+=1
                data=atoms(overlap,src,4);inputs.append(data)
                donors.extend([src]*len(data));receivers.extend([dest]*len(data))
                xy,dv=quadrature(overlap,8);phi8+=rho[src]*(dv@potential(xy))
        assert np.max(abs(coverold/vol-1))<1e-8
        assert np.max(abs(covernew/newvol-1))<1e-8
        data=np.concatenate(inputs);donors=np.array(donors);receivers=np.array(receivers)
        assert np.min(data[:,0])>0
        assert cpu.wm_set_threads(1)==1
        one=call(cpu.wm_batch,data)
        assert cpu.wm_set_threads(4)==4
        four=call(cpu.wm_batch,data);device=call(gpu.wm_cuda_batch,data)
        assert np.array_equal(one.view('u8'),four.view('u8'))
        assert np.array_equal(one.view('u8'),device.view('u8'))
        pre=np.zeros_like(ref);post=np.zeros_like(ref)
        np.add.at(pre,donors,one);np.add.at(post,receivers,one)
        polynomial_error=float(np.max(abs(pre[:,:16]-ref[:,:16])/(1+abs(ref[:,:16]))))
        assert polynomial_error<1e-12
        total_error=float(np.max(abs(post.sum(axis=0)-pre.sum(axis=0))))
        assert total_error<1e-11
        m=post[:,0];p=post[:,6:8]
        kinetic=.5*(post[:,12]+post[:,14]);heat=post[:,15]-kinetic
        pnew=(gamma-1)*heat/newvol
        assert np.all(heat>0) and np.all(m>0)
        p_error=float(np.max(abs(pnew/pressure-1)))
        assert p_error<1e-8
        mean_kinetic=.5*np.sum(p*p,axis=1)/m
        variance=kinetic-mean_kinetic
        assert np.min(variance)>-1e-14
        # Necessary covariance PSD check; explicit positive overlap atoms supply
        # a stronger geometric witness for this one remap (not a closure law).
        joint=np.zeros((len(ids),4,4));joint[:,0,0]=post[:,3];joint[:,0,1]=joint[:,1,0]=post[:,4];joint[:,1,1]=post[:,5]
        joint[:,0,2]=joint[:,2,0]=post[:,8];joint[:,0,3]=joint[:,3,0]=post[:,9]
        joint[:,1,2]=joint[:,2,1]=post[:,10];joint[:,1,3]=joint[:,3,1]=post[:,11]
        joint[:,2,2]=post[:,12];joint[:,2,3]=joint[:,3,2]=post[:,13];joint[:,3,3]=post[:,14]
        mean=np.column_stack((post[:,1:3],post[:,6:8]))/m[:,None]
        covariance=joint/m[:,None,None]-mean[:,:,None]*mean[:,None,:]
        min_eig=float(np.linalg.eigvalsh(covariance).min());assert min_eig>-1e-12
        actual_l_before=float(np.sum(pre[:,9]-pre[:,10]));actual_l_after=float(np.sum(post[:,9]-post[:,10]))
        point_before=float(np.sum(x[ids,0]*pre[:,7]-x[ids,1]*pre[:,6]))
        point_after=float(np.sum(x[ids,0]*post[:,7]-x[ids,1]*post[:,6]))
        row={'weight_fraction':fraction,'cells':len(ids),'overlap_fragments':pieces,'quadrature_atoms':len(data),
             'old_coverage_max_relative_error':float(np.max(abs(coverold/vol-1))),
             'new_coverage_max_relative_error':float(np.max(abs(covernew/newvol-1))),
             'independent_old_polynomial_moment_max_scaled_error':polynomial_error,
             'global_all_17_moment_max_absolute_change':total_error,
             'openmp_1_4_and_cuda_bitwise_equal':True,'integrated_angular_change':actual_l_after-actual_l_before,
             'point_angular_change':point_after-point_before,
             'internal_velocity_energy_over_initial_thermal':float(variance.sum()/(pressure*vol.sum()/(gamma-1))),
             'uniform_pressure_max_relative_error':p_error,'minimum_mass':float(m.min()),
             'minimum_thermal_energy':float(heat.min()),'joint_covariance_minimum_eigenvalue':min_eig,
             'potential_partition_quadrature_difference_order4':float(pre[:,16].sum()-ref[:,16].sum()),
             'potential_partition_quadrature_difference_order8':float(phi8-ref8[:,16].sum()),
             'potential_shared_atom_change':float(post[:,16].sum()-pre[:,16].sum()),
             'mass_change_fraction_max':float(np.max(abs(m/mass-1)))}
        if fraction==0: assert np.max(abs(post-ref)/(1+abs(ref)))<1e-12
        reports.append(row);print(json.dumps(row),flush=True)
        np.savez(root/f'moments_fraction{fraction:g}.npz',ids=ids,old=pre,new=post,old_volume=vol,new_volume=newvol)
    point_l=float(np.sum(x[ids,0]*mass*vel[ids,1]-x[ids,1]*mass*vel[ids,0]))
    cell_l=float(np.sum(ref[:,9]-ref[:,10]))
    point_phi=float(mass@potential(x[ids]));cell_phi=float(ref8[:,16].sum())
    result={'scope':'One finite fixed-position weight remap on canonical C checkpoint geometry. Donor fields are cell-constant with original masses, sampled circular velocities and common pressure. No hydrodynamic evolution or MPI test.',
            'state_layout':['m','qx','qy','Qxx','Qxy','Qyy','px','py','Dxx','Dxy','Dyx','Dyy','Pixx','Pixy','Piyy','E_hydro','E_potential'],
            'reports':reports,'restart_initialization':{'automatic_conversion_allowed':False,'point_L':point_l,'cell_integrated_L':cell_l,'angular_difference':cell_l-point_l,'point_potential':point_phi,'cell_integrated_potential':cell_phi,'potential_difference':cell_phi-point_phi,'mass_momentum_hydroenergy':'Preserved by constant donor field construction; geometric angular momentum and potential generally differ from old point convention.'},
            'limitations':['Positive atoms supply a realizable state for this remap only. A compact moment reconstruction for subsequent hydro/weight steps is still required.',
                           'Unresolved kinetic energy must not be hidden permanently to suppress physical shock heating; production closure and dissipation are unimplemented.',
                           'A conserved stored potential integral is not a gravity force/work evolution scheme. Independent quadrature comparisons are reported.',
                           'Zero weight change is identity in this component. Full approved fixed-mass Voronoi reduction is not implemented.',
                           'Old point-particle restart requires an explicit conservative initialization contract; this test does not authorize reinterpretation.'],
            'seconds':time.monotonic()-start,'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,code,root/'weight_remap_moments.h',root/'probe.c',root/'probe.cu',root/'cpu.so',root/'cuda.so',Path(__file__))}}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'seconds':result['seconds'],'restart_initialization':result['restart_initialization']}),flush=True)


if __name__=='__main__': main()

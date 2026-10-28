"""Actual t2 restart -> positive moment lift -> finite weight remap.

No checkpoint is rewritten and no hydrodynamic time step is taken. Geometry
clipping is the diagnostic Python reference; density/moments use shared C/CUDA.
"""
import ctypes
import hashlib
import importlib.util
import json
from pathlib import Path
import time
import numpy as np


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod


def main():
    root=Path(__file__).parent;assert not (root/'result.json').exists()
    base=Path('/gpfs/kjhan/LagForce');started=time.monotonic()
    mod=load('moments',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
    geo=load('geometry',base/'Codex_weighted_compression_geometry_20261017/audit.py')
    reader=load('reader',base/'Codex_paired_response_20261017/run.py')
    source=base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat'
    layout=json.loads((base/'Codex_paired_response_20261017/layout.json').read_text())
    raw,a,sort,time0=reader.read(source,layout);a=a[sort]
    x=np.column_stack((a['x']-3,a['y']-3));vel=np.column_stack((a['vx'],a['vy']))
    ids=np.flatnonzero(np.linalg.norm(x,axis=1)<.7);mass=a['mass'][ids];thermal=a['ie'][ids]
    assert np.all(a['w2']==0) and np.all(thermal>0)
    old=mod.polygons(x,np.zeros(len(x)),ids,geo.clip)
    volume=np.array([mod.area(p) for p in old]);rho=mass/volume;udensity=thermal/volume
    assert np.max(abs(volume/a['volume'][ids]-1))<1e-7
    target=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights']
    # Scale only conditions the four-constraint solve; coefficients are refit
    # to the ACTUAL stored velocities, thermal energies, masses and positions.
    scale=np.load(base/'Codex_weight_remap_initialization_check_v3_20261018/l2_candidate.npz')['feature_scale']
    phi0=mod.potential(x[ids]);means=[];G=np.zeros((4,4));b=np.zeros(4)
    def features(xy,k):
        r=xy-x[ids[k]];v=vel[ids[k]]
        return np.column_stack((r,r[:,0]*v[1]-r[:,1]*v[0],mod.potential(xy)-phi0[k]))/scale
    for k,p in enumerate(old):
        xy,dv=mod.quadrature(p,8);f=features(xy,k);mean=dv@f/volume[k]
        g=rho[k]*(f-mean);G+=g.T@(dv[:,None]*g);b-=rho[k]*(dv@f);means.append(mean)
    lam=np.linalg.solve(G,b);means=np.array(means);coeff=lam/scale
    cpu=ctypes.CDLL(str(root/'cpu.so'));gpu=ctypes.CDLL(str(root/'cuda.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for fun in (cpu.wm_lift_batch,gpu.wm_lift_cuda_batch):
        fun.argtypes=[ctypes.c_size_t,arr,arr];fun.restype=ctypes.c_int
    def call(fun,data):
        data=np.ascontiguousarray(data);out=np.empty((len(data),17))
        assert fun(len(data),data,out)==0;assert np.isfinite(out).all();return out
    def inputs(poly,k,order):
        xy,dv=mod.quadrature(poly,order);n=len(xy)
        return np.column_stack((dv,xy,np.tile(x[ids[k]],(n,1)),np.tile(vel[ids[k]],(n,1)),
                                np.full(n,rho[k]),np.full(n,udensity[k]),mod.potential(xy),
                                np.full(n,phi0[k]),np.tile(coeff,(n,1)),np.full(n,means[k]@lam)))
    assert cpu.wm_set_threads(1)==1
    reference={}
    for order in (4,8):
        reference[order]=np.array([call(cpu.wm_lift_batch,inputs(p,k,order)).sum(axis=0) for k,p in enumerate(old)])
    point=np.array([mass.sum(),np.sum(mass*vel[ids,0]),np.sum(mass*vel[ids,1]),
                    np.sum(mass*(x[ids,0]*vel[ids,1]-x[ids,1]*vel[ids,0])),
                    np.sum(thermal+.5*mass*np.sum(vel[ids]**2,axis=1)),mass@phi0])
    def invariants(q):
        return np.array([q[:,0].sum(),q[:,6].sum(),q[:,7].sum(),np.sum(q[:,9]-q[:,10]),q[:,15].sum(),q[:,16].sum()])
    lo=np.array([p.min(axis=0) for p in old]);hi=np.array([p.max(axis=0) for p in old])
    reports=[]
    for fraction in (.25,1.):
        new=mod.polygons(x,fraction*target,ids,geo.clip);newvolume=np.array([mod.area(p) for p in new])
        fragments=[];coverold=np.zeros(len(ids));covernew=np.zeros(len(ids))
        for dest,p in enumerate(new):
            candidates=np.flatnonzero(np.all(hi>=p.min(axis=0),axis=1)&np.all(lo<=p.max(axis=0),axis=1))
            for src in candidates:
                overlap=old[src].copy()
                for aa,bb in zip(p,np.roll(p,-1,axis=0)):
                    edge=bb-aa;normal=np.array([edge[1],-edge[0]])
                    overlap=geo.clip(overlap,normal,normal@aa)
                    if len(overlap)<3:break
                av=mod.area(overlap)
                if av<=1e-18:continue
                fragments.append((src,dest,overlap));coverold[src]+=av;covernew[dest]+=av
        assert np.max(abs(coverold/volume-1))<1e-8 and np.max(abs(covernew/newvolume-1))<1e-8
        for order in (4,8):
            pre=np.zeros((len(ids),17));post=np.zeros_like(pre);expected_thermal=np.zeros(len(ids))
            min_factor=np.inf;max_oracle=0.;atoms=0
            # Batch chunks keep the GPU footprint small; OMP and GPU execute
            # the same density evaluation and 17-moment contribution.
            for begin in range(0,len(fragments),128):
                batch=[];donors=[];receivers=[]
                for src,dest,poly in fragments[begin:begin+128]:
                    data=inputs(poly,src,order);batch.append(data)
                    donors.extend([src]*len(data));receivers.extend([dest]*len(data))
                    expected_thermal[dest]+=udensity[src]*mod.area(poly)
                data=np.concatenate(batch);donors=np.array(donors);receivers=np.array(receivers)
                assert cpu.wm_set_threads(1)==1;one=call(cpu.wm_lift_batch,data)
                assert cpu.wm_set_threads(4)==4;four=call(cpu.wm_lift_batch,data)
                device=call(gpu.wm_lift_cuda_batch,data)
                assert np.array_equal(one.view('u8'),four.view('u8'))
                assert np.array_equal(one.view('u8'),device.view('u8'))
                rel=data[:,1:3]-data[:,3:5]
                feature=(rel[:,0]*coeff[0]+rel[:,1]*coeff[1]+(rel[:,0]*data[:,6]-rel[:,1]*data[:,5])*coeff[2]
                         +(data[:,9]-data[:,10])*coeff[3]-data[:,15])
                factor=1+data[:,7]*feature;oracle=data[:,0]*data[:,7]*factor
                max_oracle=max(max_oracle,float(np.max(abs(one[:,0]-oracle))))
                assert max_oracle<1e-14
                min_factor=min(min_factor,float(factor.min()));atoms+=len(data)
                np.add.at(pre,donors,one);np.add.at(post,receivers,one)
            heat=post[:,15]-.5*(post[:,12]+post[:,14]);assert np.all(heat>0)
            conservation=invariants(post)-invariants(pre)
            assert np.max(abs(conservation))<1e-11
            p_error=float(np.max(abs(heat/expected_thermal-1)));assert p_error<1e-8
            var=.5*(post[:,12]+post[:,14])-.5*np.sum(post[:,6:8]**2,axis=1)/post[:,0]
            assert np.min(var)>-1e-14
            row={'weight_fraction':fraction,'order':order,'fragments':len(fragments),'atoms':atoms,
                 'CPU_OMP1_OMP4_CUDA_bitwise_equal':True,'mass_oracle_max_absolute_error':max_oracle,
                 'minimum_sampled_density_factor':min_factor,'all_mass_and_thermal_positive':True,
                 'shared_atom_invariant_changes':conservation.tolist(),
                 'independent_donor_quadrature_invariant_difference':(invariants(pre)-invariants(reference[order])).tolist(),
                 'post_minus_original_point_invariants':(invariants(post)-point).tolist(),
                 'thermal_overlap_max_relative_error':p_error,
                 'unresolved_kinetic_over_initial_thermal':float(var.sum()/thermal.sum())}
            reports.append(row);print(json.dumps(row),flush=True)
            np.savez(root/f'state_w{fraction:g}_q{order}.npz',ids=ids,state=post,volume=newvolume)
    # Nonpositive state must be rejected, not floored.
    bad=inputs(old[0],0,4)[:1].copy();bad[:,8]=-1;out=np.empty((1,17))
    assert cpu.wm_lift_batch(1,bad,out)!=0 and gpu.wm_lift_cuda_batch(1,bad,out)!=0
    result={'scope':'Actual t2 restart fields, central r<.7 patch of1590cells only. Conditional initialization plus finite weight remapping; no physical time evolution or MPI.',
            'source':str(source),'source_time':time0,'invariant_layout':['M','Px','Py','Lz','E_hydro','E_potential'],
            'initialization_error_order8':(invariants(reference[8])-point).tolist(),
            'original_thermal_min':float(thermal.min()),'original_thermal_sum':float(thermal.sum()),
            'original_pressure_range':[float(a['pressure'][ids].min()),float(a['pressure'][ids].max())],
            'reports':reports,'negative_thermal_input_rejected_cpu_cuda':True,
            'remaining':['Continuous positivity for these newly fitted actual-velocity coefficients remains to be bounded; positive quadrature atoms alone do not prove it.',
                         'Retaining moments and potential does not define a pressure/gravity evolution closure or generator motion.',
                         'No production checkpoint changed; approved zero-weight fixed-mass reduction remains a full-solver requirement.'],
            'seconds':time.monotonic()-started,
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,root/'weight_remap_moments.h',root/'probe.c',root/'probe.cu',root/'cpu.so',root/'cuda.so',Path(__file__))}}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(root/'lift_coefficients.npz',ids=ids,coefficient=lam,feature_scale=scale,cell_feature_mean=means,volume=volume)
    print(json.dumps({'seconds':result['seconds'],'initialization':result['initialization_error_order8']}),flush=True)


if __name__=='__main__':main()

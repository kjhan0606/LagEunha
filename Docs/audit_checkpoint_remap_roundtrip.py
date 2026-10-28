"""Actual checkpoint fragment remap 0 -> full weights -> quarter -> 0.

Keep donor density/velocity/thermal fields on every retained fragment. This is
a finite pure-remapping audit, not a compact-state hydrodynamic integrator.
"""
import ctypes
import hashlib
import importlib.util
import json
from pathlib import Path
import time
import numpy as np


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def main():
    root=Path(__file__).parent;assert not (root/'result.json').exists()
    base=Path('/gpfs/kjhan/LagForce');kernel=base/'Codex_checkpoint_weight_remap_20261018'
    mod=load('moments',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
    geo=load('geometry',base/'Codex_weighted_compression_geometry_20261017/audit.py')
    reader=load('reader',base/'Codex_paired_response_20261017/run.py')
    source=base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat'
    layout=json.loads((base/'Codex_paired_response_20261017/layout.json').read_text())
    _,a,sort,t=reader.read(source,layout);a=a[sort]
    x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']))
    z=np.load(kernel/'lift_coefficients.npz');ids=z['ids'];vol=z['volume']
    mass=a['mass'][ids];thermal=a['ie'][ids];rho=mass/vol
    coeff=z['coefficient']/z['feature_scale'];mean=z['cell_feature_mean']@z['coefficient']
    phi0=mod.potential(x[ids]);target=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights']
    old=mod.polygons(x,np.zeros(len(x)),ids,geo.clip)
    cpu=ctypes.CDLL(str(kernel/'cpu.so'));gpu=ctypes.CDLL(str(kernel/'cuda.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for fun in (cpu.wm_lift_batch,gpu.wm_lift_cuda_batch):
        fun.argtypes=[ctypes.c_size_t,arr,arr];fun.restype=ctypes.c_int
    assert cpu.wm_set_threads(4)==4
    def inputs(poly,donor):
        xy,dv=mod.quadrature(poly,8);n=len(xy);k=donor
        return np.column_stack((dv,xy,np.tile(x[ids[k]],(n,1)),np.tile(v[ids[k]],(n,1)),
                                np.full(n,rho[k]),np.full(n,thermal[k]/vol[k]),mod.potential(xy),
                                np.full(n,phi0[k]),np.tile(coeff,(n,1)),np.full(n,mean[k])))
    def integrate(fragments):
        result=np.zeros((len(ids),17));minimum=np.inf;atoms=0
        for begin in range(0,len(fragments),128):
            blocks=[];owner=[]
            for donor,cell,p in fragments[begin:begin+128]:
                data=inputs(p,donor);blocks.append(data);owner.extend([cell]*len(data))
            data=np.ascontiguousarray(np.concatenate(blocks));one=np.empty((len(data),17));device=np.empty_like(one)
            assert cpu.wm_lift_batch(len(data),data,one)==0
            assert gpu.wm_lift_cuda_batch(len(data),data,device)==0
            assert np.array_equal(one.view('u8'),device.view('u8'))
            assert np.isfinite(one).all() and np.all(one[:,0]>0)
            minimum=min(minimum,float(np.min(one[:,0]/(data[:,0]*data[:,7]))))
            np.add.at(result,np.array(owner),one);atoms+=len(data)
        return result,minimum,atoms
    def invariants(q):
        return np.array([q[:,0].sum(),q[:,6].sum(),q[:,7].sum(),np.sum(q[:,9]-q[:,10]),q[:,15].sum(),q[:,16].sum()])
    fragments=[(k,k,p) for k,p in enumerate(old)];before,_,_=integrate(fragments);start=time.monotonic()
    rows=[]
    for fraction in (1.,.25,0.):
        new=mod.polygons(x,fraction*target,ids,geo.clip)
        lo=np.array([p.min(axis=0) for p in new]);hi=np.array([p.max(axis=0) for p in new])
        newvol=np.array([mod.area(p) for p in new]);coverage=np.zeros(len(ids));out=[];fragment_error=0.;discard=0.
        for donor,owner,p in fragments:
            area=mod.area(p);covered=0.
            candidates=np.flatnonzero(np.all(hi>=p.min(axis=0),axis=1)&np.all(lo<=p.max(axis=0),axis=1))
            for dest in candidates:
                overlap=p.copy();receiver=new[dest]
                for aa,bb in zip(receiver,np.roll(receiver,-1,axis=0)):
                    edge=bb-aa;normal=np.array([edge[1],-edge[0]])
                    overlap=geo.clip(overlap,normal,normal@aa)
                    if len(overlap)<3:break
                part=mod.area(overlap)
                if part<=1e-18:
                    discard+=max(part,0.);continue
                coverage[dest]+=part;covered+=part;out.append((donor,dest,overlap))
            fragment_error=max(fragment_error,abs(covered-area))
        assert np.max(abs(coverage/newvol-1))<1e-8
        fragments=out;state,minfactor,atoms=integrate(fragments)
        heat=state[:,15]-.5*(state[:,12]+state[:,14]);assert np.all(heat>0)
        invdiff=invariants(state)-invariants(before);assert np.max(abs(invdiff))<1e-10
        variance=.5*(state[:,12]+state[:,14])-.5*np.sum(state[:,6:8]**2,axis=1)/state[:,0]
        row={'weight_fraction':fraction,'fragments':len(fragments),'quadrature_atoms':atoms,
             'CPU_OMP4_CUDA_bitwise_equal':True,'minimum_sampled_density_factor':minfactor,
             'maximum_fragment_area_closure_error':fragment_error,
             'discarded_positive_sliver_area':discard,'receiver_coverage_max_relative_error':float(np.max(abs(coverage/newvol-1))),
             'aggregate_invariant_change_from_start':invdiff.tolist(),
             'unresolved_kinetic_over_initial_thermal':float(variance.sum()/thermal.sum())}
        if fraction==0:
            row.update(return_mass_max_relative_error=float(np.max(abs(state[:,0]/mass-1))),
                       return_thermal_max_relative_error=float(np.max(abs(heat/thermal-1))),
                       return_velocity_max_absolute_error=float(np.max(abs(state[:,6:8]/state[:,0,None]-v[ids]))),
                       return_all17_moments_max_absolute_error=float(np.max(abs(state-before))))
            assert row['return_mass_max_relative_error']<1e-8
            assert row['return_thermal_max_relative_error']<1e-8
            assert row['return_velocity_max_absolute_error']<1e-10
        rows.append(row);print(json.dumps(row),flush=True)
    np.savez(root/'final.npz',ids=ids,initial=before,final=state)
    result={'scope':'Actual t2 central1590cell patch, fixed generators. Sequential retained-fragment remap 0 ->1 ->.25 ->0. No physical time, pressure/gravity update, compact closure or MPI.',
            'source':str(source),'source_time':t,'invariant_layout':['M','Px','Py','Lz','E_hydro','E_potential'],
            'rows':rows,'seconds':time.monotonic()-start,
            'zero_weight_interpretation':'Return uses intersections and original donor fields carried through prior fragments. No mass reset, thermal correction or discarded velocity moments. Round-trip identity here is a pure-remap component test, not the complete fixed-mass Voronoi evolution limit.',
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,kernel/'lift_coefficients.npz',kernel/'weight_remap_moments.h',kernel/'cpu.so',kernel/'cuda.so',Path(__file__))}}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()

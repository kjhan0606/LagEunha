"""Actual material moment drift versus support in newly moved power cells.

This intentionally omits moving-face moment flux, to determine when it cannot
possibly be repaired by a positive density reconstruction alone. Negative
expectation of a product of two nonnegative cell halfspace distances is a
necessary-support failure, independent of polynomial density degree.
No pressure/gravity evolution, remap or physical weight law is claimed.
"""
import ctypes
import hashlib
import importlib.util
import json
from pathlib import Path
import time
import numpy as np

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def support(poly,q):
    # Work in normalized centered coordinates to reduce cancellation.
    origin=poly.mean(axis=0);p=poly-origin;edge=np.roll(p,-1,axis=0)-p
    area=.5*np.sum(p[:,0]*np.roll(p,-1,axis=0)[:,1]-p[:,1]*np.roll(p,-1,axis=0)[:,0]);h=np.sqrt(area);p/=h
    edge=np.roll(p,-1,axis=0)-p;n=np.column_stack((edge[:,1],-edge[:,0]));n/=np.linalg.norm(n,axis=1)[:,None]
    b=np.sum(n*p,axis=1);mean=q[1:3]/q[0];c=(mean-origin)/h
    Q=np.array([[q[3],q[4]],[q[4],q[5]]])/q[0];Cov=(Q-np.outer(mean,mean))/(h*h)
    gap=b-n@c
    product=np.outer(gap,gap)+n@Cov@n.T
    return float(gap.min()),float(product.min())

def main():
    root=Path(__file__).parent;assert not (root/'result.json').exists();started=time.monotonic();base=Path('/gpfs/kjhan/LagForce')
    modpath=base/'Codex_weight_remap_moments_centered_20261017/audit.py';geopath=base/'Codex_weighted_compression_geometry_20261017/audit.py';readpath=base/'Codex_paired_response_20261017/run.py';decodepath=base/'Codex_moment_pressure_kick_20261018/audit.py'
    mod=load('moments',modpath);geo=load('geo',geopath);reader=load('reader',readpath);decoder=load('decoder',decodepath)
    source=base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat';layout=base/'Codex_paired_response_20261017/layout.json'
    _,a,sort,t=reader.read(source,json.loads(layout.read_text()));a=a[sort]
    x=np.column_stack((a['x']-3,a['y']-3));vel=np.column_stack((a['vx'],a['vy']));targetpath=base/'Codex_evolving_weight_transport_20261017/weights_0.npz';target=np.load(targetpath)['weights']
    cpu=ctypes.CDLL(str(root/'cpu.so'));gpu=ctypes.CDLL(str(root/'cuda.so'));arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for f in (cpu.md_batch,gpu.md_cuda_batch):f.argtypes=[ctypes.c_size_t,arr,arr];f.restype=ctypes.c_int
    def call(f,q,dt):
        data=np.ascontiguousarray(np.column_stack((q[:,:16],np.full(len(q),dt))));out=np.empty((len(q),16));assert f(len(q),data,out)==0;return out
    rows=[];paths=[Path(__file__),modpath,geopath,readpath,decodepath,source,layout,targetpath,root/'moment_material_drift.h',root/'probe.c',root/'probe.cu',root/'cpu.so',root/'cuda.so']
    for fraction in (.25,1.):
        statepath=base/f'Codex_checkpoint_weight_remap_20261018/state_w{fraction:g}_q8.npz';paths.append(statepath);saved=np.load(statepath);ids=saved['ids'];q=saved['state'];V=saved['volume']
        c,v,S,C,B,R,U=decoder.decode(q);generator=vel.copy();generator[ids]=v+np.einsum('nij,nj->ni',B,x[ids]-c)
        for dt in (0.,.001,.005,.02,.05,.1):
            cpu.md_threads(1);one=call(cpu.md_batch,q,dt);cpu.md_threads(4);four=call(cpu.md_batch,q,dt);device=call(gpu.md_cuda_batch,q,dt)
            assert np.array_equal(one.view('u8'),four.view('u8')) and np.array_equal(one.view('u8'),device.view('u8'))
            back=call(cpu.md_batch,one,-dt);inverse_error=float(np.max(abs(back-q[:,:16])));assert inverse_error<1e-12
            angular=float(np.sum((one[:,9]-one[:,10])-(q[:,9]-q[:,10])));assert abs(angular)<1e-12
            assert np.array_equal(one[:,[0,6,7,12,13,14,15]],q[:,[0,6,7,12,13,14,15]])
            row={'weight_fraction':fraction,'dt':dt,'CPU_OMP1_OMP4_CUDA_bitwise_equal':True,'inverse_drift_max_absolute_error':inverse_error,'angular_change':angular,
                 'mass_momentum_kinetic_hydro_energy_unchanged':True,'omitted_R_position_covariance_trace_over_initial_trace_max':float(np.max(dt*dt*np.trace(R,axis1=1,axis2=2)/np.trace(S,axis1=1,axis2=2)))}
            try:polys=mod.polygons(x+dt*generator,fraction*target,ids,geo.clip)
            except (AssertionError,RuntimeError) as error:
                row.update(geometry_valid=False,geometry_failure=str(error));rows.append(row);print(json.dumps(row),flush=True);continue
            bounds=np.array([support(p,z) for p,z in zip(polys,one)]);bad=np.any(bounds < -1e-9,axis=1)
            row.update(geometry_valid=True,minimum_normalized_mean_gap=float(bounds[:,0].min()),minimum_normalized_halfspace_product=float(bounds[:,1].min()),
                       necessary_support_failures=int(bad.sum()),failed_ids=ids[bad].tolist())
            if dt==0:assert not bad.any()
            rows.append(row);print(json.dumps(row),flush=True)
            np.savez(root/f'drift_w{fraction:g}_dt{dt:g}.npz',ids=ids,hydro_moments=one,sites=x+dt*generator,support_bounds=bounds)
    result={'scope':'Full material free drift of actual finite-remap moments; generator moves with affine reconstructed point velocity in central patch, prescribed saved velocity outside. Fixed weight differences; no force, remap or moment face flux.',
            'source_time':t,'rows':rows,'support_criterion':'Within a convex polygon every inward halfspace distance l_i>=0, hence E[l_i]>=0 and E[l_i*l_j]>=0. Values below -1e-9 in normalized coordinates are recorded as necessary-condition failures. Passing is not sufficient for density realizability.',
            'gravity':'Potential is excluded from the16-output drift; unchanged potential would not be physically valid for moving mass in external gravity.',
            'lineage':'2026-10-09 positive material-fragment drift existed; this is its compact C/OpenMP/CUDA full-moment counterpart and actual canonical geometry support test.',
            'limits':['A material packet is not automatically the contents of a moving power cell. Missing face moments must not be hidden by density clipping or higher polynomial degree.', 'Affine generator velocity omits a spatial realization of residual motion. No general closure or zero-weight reduction proof.', 'No long-time Kepler or MPI validation in this component.'],
            'seconds':time.monotonic()-started,'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()

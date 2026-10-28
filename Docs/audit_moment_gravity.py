"""Shared C/CUDA material moment source and independent spline-work checks."""
import ctypes
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def main():
    root=Path(__file__).parent;base=Path('/gpfs/kjhan/LagForce')
    assert not (root/'result.json').exists()
    mod=load('moments',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
    geo=load('geometry',base/'Codex_weighted_compression_geometry_20261017/audit.py')
    reader=load('reader',base/'Codex_paired_response_20261017/run.py')
    source=base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat'
    layout=json.loads((base/'Codex_paired_response_20261017/layout.json').read_text())
    _,a,sort,t=reader.read(source,layout);a=a[sort]
    x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']))
    fit=base/'Codex_checkpoint_weight_remap_20261018/lift_coefficients.npz'
    z=np.load(fit);ids=z['ids'];volume=z['volume'];rho=a['mass'][ids]/volume
    lam=z['coefficient'];scale=z['feature_scale'];means=z['cell_feature_mean']
    phi0=mod.potential(x[ids]);polys=mod.polygons(x,np.zeros(len(x)),ids,geo.clip)
    batches=[];owner=[]
    for k,p in enumerate(polys):
        xy,dv=mod.quadrature(p,8);r=xy-x[ids[k]];vel=v[ids[k]]
        f=np.column_stack((r,r[:,0]*vel[1]-r[:,1]*vel[0],mod.potential(xy)-phi0[k]))/scale
        density=rho[k]*(1+rho[k]*((f-means[k])@lam));assert np.all(density>0)
        batches.append(np.column_stack((density*dv,xy,np.tile(vel,(len(xy),1)),np.full(len(xy),.5))))
        owner.extend([k]*len(xy))
    data=np.ascontiguousarray(np.concatenate(batches));owner=np.array(owner)
    cpu=ctypes.CDLL(str(root/'cpu.so'));gpu=ctypes.CDLL(str(root/'cuda.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for fun in (cpu.wg_batch,gpu.wg_cuda_batch):
        fun.argtypes=[ctypes.c_size_t,arr,arr];fun.restype=ctypes.c_int
    def call(fun,d):
        d=np.ascontiguousarray(d);out=np.empty((len(d),20));assert fun(len(d),d,out)==0
        assert np.isfinite(out).all();return out
    assert cpu.wg_threads(1)==1;one=call(cpu.wg_batch,data)
    assert cpu.wg_threads(4)==4;four=call(cpu.wg_batch,data);device=call(gpu.wg_cuda_batch,data)
    assert np.array_equal(one.view('u8'),four.view('u8'))
    parity=float(np.max(abs(device-one)/(1+abs(one))));assert parity<1e-13
    assert np.all(one[:,15]+one[:,16]==0)
    phi_error=float(np.max(abs(one[:,19]-mod.potential(data[:,1:3]))));assert phi_error<1e-13
    thermal_rate=(one[:,15]-.5*(one[:,12]+one[:,14]))/data[:,0]
    torque_rate=(one[:,9]-one[:,10])/data[:,0]
    assert np.max(abs(thermal_rate))<1e-13 and np.max(abs(torque_rate))<1e-13
    # Include r=0, both spline joins and outer points in a separate device check.
    radius=np.array([0.,1e-8,.2499999,.25,.2500001,.4999999,.5,.5000001,1.,3.])
    boundary=np.column_stack((np.ones(len(radius)),radius,np.zeros(len(radius)),
                              np.full(len(radius),.3),np.full(len(radius),-.7),np.full(len(radius),.5)))
    bd=call(gpu.wg_cuda_batch,boundary)
    assert np.max(abs(bd[:,19]-mod.potential(boundary[:,1:3])))<1e-13
    assert np.all(bd[0,17:19]==0)

    # Independent centered difference of the moment definition along material
    # tangent xdot=v, vdot=g; NOT a physical time integrator or fluid solution.
    sample=data[::max(1,len(data)//2048)].copy();sample[:,0]=1
    exact=call(cpu.wg_batch,sample)[:,:17]
    gravity=call(cpu.wg_batch,sample)[:,17:19]
    def values(pos,vel):
        xx,yy=pos.T;vx,vy=vel.T
        return np.column_stack((np.ones(len(pos)),xx,yy,xx*xx,xx*yy,yy*yy,vx,vy,
                                xx*vx,xx*vy,yy*vx,yy*vy,vx*vx,vx*vy,vy*vy,
                                .5*(vx*vx+vy*vy),mod.potential(pos)))
    differences=[]
    for eps in (1e-4,5e-5,2.5e-5):
        plus=values(sample[:,1:3]+eps*sample[:,3:5],sample[:,3:5]+eps*gravity)
        minus=values(sample[:,1:3]-eps*sample[:,3:5],sample[:,3:5]-eps*gravity)
        fd=(plus-minus)/(2*eps)
        differences.append({'epsilon':eps,'first16_max_scaled_error':float(np.max(abs(fd[:,:16]-exact[:,:16])/(1+abs(exact[:,:16])))),
                            'potential_rate_max_scaled_error':float(np.max(abs(fd[:,16]-exact[:,16])/(1+abs(exact[:,16]))))})
    assert differences[-1]['first16_max_scaled_error']<1e-9
    assert differences[-1]['potential_rate_max_scaled_error']<1e-6
    cell=np.zeros((len(ids),17));np.add.at(cell,owner,one[:,:17])
    pointdata=np.column_stack((a['mass'][ids],x[ids],v[ids],np.full(len(ids),.5)))
    point=call(cpu.wg_batch,pointdata)
    da=(cell[:,6:8]-point[:,6:8])/a['mass'][ids,None]
    difference=float(np.sqrt(np.average(np.sum(da*da,axis=1),weights=a['mass'][ids])))
    denominator=float(np.sqrt(np.average(np.sum(point[:,17:19]**2,axis=1),weights=a['mass'][ids])))
    result={'scope':'Actual t2 central1590cell positive density representation, zero weights. Material-quadrature gravity/kinematic source only; no pressure, remap, boundary flux, finite step or MPI.',
            'source':str(source),'source_time':t,'atoms':len(data),'cpu_omp1_4_bitwise_equal':True,
            'cpu_cuda_bitwise_equal':bool(np.array_equal(one.view('u8'),device.view('u8'))),
            'cpu_cuda_max_scaled_difference':parity,'spline_potential_numpy_max_absolute_difference':phi_error,
            'maximum_thermal_source_per_mass_residual':float(np.max(abs(thermal_rate))),
            'maximum_central_torque_per_mass_residual':float(np.max(abs(torque_rate))),
            'total_hydro_plus_potential_source':float(np.sum(cell[:,15]+cell[:,16])),
            'finite_difference_checks':differences,
            'zero_weight_force_comparison':{'integrated_vs_point_acceleration_mass_weighted_RMS':difference,
                                          'fraction_of_point_acceleration_RMS':difference/denominator,
                                          'integrated_gravity_power':float(cell[:,15].sum()),'point_gravity_power':float(point[:,15].sum())},
            'interpretation':'Correct moment gravity work does not imply exact recovery of the old point-force baseline. The measured zero-weight force difference must be addressed by the full reduction contract, not silently enabled.',
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,fit,root/'moment_spline_gravity.h',root/'probe.c',root/'probe.cu',root/'cpu.so',root/'cuda.so',Path(__file__))}}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()

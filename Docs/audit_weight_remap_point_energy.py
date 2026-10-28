"""Actual t2 point-state budget, no subcell gravity quadrature or evolution.

Independent variable s is weight amplitude: w(s)=s*w_target, evaluated at s=0.
All reported rates are per unit s, NOT per physical time. Geometry is a Python
reference; the pair kernel runs on C/OpenMP and actual CUDA without fallback.
"""
import ctypes
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def main():
    root=Path(__file__).parent
    assert not (root/'result.json').exists()
    base=Path('/gpfs/kjhan/LagForce')
    reader_path=base/'Codex_paired_response_20261017/run.py'
    geo_path=base/'Codex_weighted_compression_geometry_20261017/audit.py'
    potential_path=base/'Codex_weight_remap_moments_centered_20261017/audit.py'
    reader=load('reader',reader_path);geo=load('geo',geo_path);pot=load('pot',potential_path)
    source=base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat'
    layout_path=base/'Codex_paired_response_20261017/layout.json'
    _,a,sort,t=reader.read(source,json.loads(layout_path.read_text()));a=a[sort]
    x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']))
    ids=np.flatnonzero(np.linalg.norm(x,axis=1)<.7)
    volume,_,lap,_=geo.geometry(x,np.zeros(len(a)),v,ids)
    volume_error=float(np.max(abs(volume/a['volume'][ids]-1)))
    assert volume_error<1e-7
    target_path=base/'Codex_evolving_weight_transport_20261017/weights_0.npz'
    target=np.load(target_path)['weights']
    ll=lap.tocoo();ii,jj=ids[ll.row],ll.col
    take=(ii<jj)&(target[ii]!=target[jj])
    ii,jj,c=ii[take],jj[take],-ll.data[take]
    assert np.isin(ii,ids).all() and np.isin(jj,ids).all()
    swept=c*(target[ii]-target[jj])
    rho=a['mass']/a['volume'];rate=swept*np.where(swept>=0,rho[jj],rho[ii])
    phi=pot.potential(x);e=a['ie']/a['mass']
    data=np.ascontiguousarray(np.column_stack((rate,v[ii],v[jj],e[ii],e[jj],phi[ii],phi[jj])))
    cpu=ctypes.CDLL(str(root/'cpu.so'));gpu=ctypes.CDLL(str(root/'cuda.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for fun in (cpu.wp_batch,gpu.wp_cuda_batch):
        fun.argtypes=[ctypes.c_size_t,arr,arr];fun.restype=ctypes.c_int
    def call(fun,d):
        out=np.empty((len(d),6));assert fun(len(d),d,out)==0
        assert np.isfinite(out).all();return out
    assert cpu.wp_threads(1)==1;one=call(cpu.wp_batch,data)
    assert cpu.wp_threads(4)==4;four=call(cpu.wp_batch,data)
    device=call(gpu.wp_cuda_batch,data)
    assert np.array_equal(one.view('u8'),four.view('u8'))
    assert np.array_equal(one.view('u8'),device.view('u8'))
    zero=data.copy();zero[:,0]=0
    assert np.all(call(cpu.wp_batch,zero)==0) and np.all(call(gpu.wp_cuda_batch,zero)==0)
    shifted=data.copy();shifted[:,7:9]+=8.
    shifted_result=call(cpu.wp_batch,shifted)
    gauge_error=float(np.max(abs(shifted_result-one)))
    assert gauge_error<1e-14
    rhs=np.zeros((len(a),4));np.add.at(rhs,ii,one[:,:4]);np.add.at(rhs,jj,-one[:,:4])
    correction=np.zeros(len(a));np.add.at(correction,ii,one[:,4]);np.add.at(correction,jj,one[:,5])
    kdot=np.sum(v*rhs[:,1:3],axis=1)-.5*np.sum(v*v,axis=1)*rhs[:,0]
    udot=rhs[:,3]-kdot
    phidot=float(phi@rhs[:,0]);corrected=udot+correction
    total_residual=float(np.sum(rhs[:,3]+correction)+phidot)
    assert abs(total_residual)<1e-13 and abs(kdot.sum())<1e-13
    U=float(a['ie'][ids].sum())
    regions={}
    r=np.linalg.norm(x,axis=1)
    for name,mask in [('inner',r<.5),('guard',(r>=.5)&(r<.7)),('patch',r<.7)]:
        rel=correction[mask]/a['ie'][mask]
        regions[name]={'source_sum_over_region_U_per_s':float(correction[mask].sum()/a['ie'][mask].sum()),
                       'source_abs_sum_over_region_U_per_s':float(abs(correction[mask]).sum()/a['ie'][mask].sum()),
                       'source_over_cell_U_quantiles_per_s':np.quantile(rel,[0,.01,.5,.99,1]).tolist()}
    pair_potential=rate*(phi[ii]-phi[jj])
    pair_energy=one[:,4]+one[:,5]+pair_potential
    assert np.max(abs(pair_energy))<1e-14
    result={'scope':'Actual t2 saved mass/vx/vy/ie at zero weights; weight-amplitude derivative only. No physical time integration, no pressure force, no MPI, no production changes.',
            'source':str(source),'source_time':t,'patch_cells':len(ids),'active_faces':len(ii),
            'parameter':'w(s)=s*w_target at s=0; rates per unit s, not per physical time',
            'CPU_OMP1_OMP4_CUDA_bitwise_equal':True,'zero_relative_weight_rate_flux_zero':True,
            'point_potential_gauge_shift_max_error':gauge_error,'C_volume_max_relative_error':volume_error,
            'patch_initial_U':U,'uncorrected_point_total_energy_rate':float(rhs[:,3].sum()+phidot),
            'corrected_point_total_energy_rate':total_residual,'kinetic_rate':float(kdot.sum()),
            'added_thermal_rate_over_patch_U_per_s':float(correction.sum()/U),
            'added_absolute_thermal_rate_over_patch_U_per_s':float(abs(correction).sum()/U),
            'point_angular_rate':float(np.sum(x[:,0]*rhs[:,2]-x[:,1]*rhs[:,1])),
            'pair_energy_residual_max':float(np.max(abs(pair_energy))),'regions':regions,
            'interpretation':'Cheap face potential correction closes the chosen point total-energy budget but creates a thermal source in pure weight relabeling. Total-energy conservation alone does not justify adopting this model. Point angular defect is unchanged. Requires a compact compatible representation/closure, not dense cell gravity quadrature.',
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,reader_path,layout_path,geo_path,potential_path,target_path,Path(__file__),root/'weight_remap_flux.h',root/'weight_remap_point_energy.h',root/'probe.c',root/'probe.cu',root/'cpu.so',root/'cuda.so']}}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()

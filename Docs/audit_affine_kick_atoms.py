"""Independent positive-atom oracle for the finite affine kick component."""
import ctypes
import hashlib
import json
from pathlib import Path
import numpy as np

def moments(m,x,v,thermal,phi):
    return np.array([m.sum(),*(m@x),m@(x[:,0]**2),m@(x[:,0]*x[:,1]),m@(x[:,1]**2),
                     *(m@v),*(x.T@(m[:,None]*v)).ravel(),m@(v[:,0]**2),m@(v[:,0]*v[:,1]),m@(v[:,1]**2),
                     thermal+.5*m@np.sum(v*v,axis=1),m@phi])

def main():
    root=Path(__file__).parent;assert not (root/'atomic_validation.json').exists()
    cpu=ctypes.CDLL(str(root/'cpu.so'));gpu=ctypes.CDLL(str(root/'cuda.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for f in (cpu.mk_batch,gpu.mk_cuda_batch):f.argtypes=[ctypes.c_size_t,arr,arr];f.restype=ctypes.c_int
    def call(f,data):
        out=np.empty((1,17));assert f(1,np.ascontiguousarray(data[None,:]),out)==0;return out[0]
    rng=np.random.default_rng(20261018);x=rng.uniform(-.3,.3,(31,2));m=rng.uniform(.1,1,len(x));v=rng.normal(size=x.shape)
    thermal=3.;phi=-1/np.sqrt(np.sum(x*x,axis=1)+.25);q=moments(m,x,v,thermal,phi)
    c=m@x/m.sum();A=np.array([[.7,-1.3],[.2,-.8]]);acc=np.array([.3,-.5])+(x-c)@A.T
    F=m@acc;J=x.T@(m[:,None]*acc);rows=[]
    for dt in (0.,.001,.1,1.,-1.):
        newv=v+dt*acc;oracle=moments(m,x,newv,thermal,phi);dE=oracle[15]-q[15]
        data=np.r_[q,F,J.ravel(),dt,dE];cpu.mk_threads(1);one=call(cpu.mk_batch,data);cpu.mk_threads(4);four=call(cpu.mk_batch,data);device=call(gpu.mk_cuda_batch,data)
        assert np.array_equal(one.view('u8'),four.view('u8')) and np.array_equal(one.view('u8'),device.view('u8'))
        error=float(np.max(abs(one-oracle)));assert error<1e-12
        reverse=call(cpu.mk_batch,np.r_[one,F,J.ravel(),-dt,-dE]);rev=float(np.max(abs(reverse-q)));assert rev<1e-12
        U=one[15]-.5*(one[12]+one[14]);assert abs(U-thermal)<1e-12
        rows.append({'dt':dt,'CPU_OMP1_OMP4_CUDA_bitwise_equal':True,'direct_atom_max_absolute_error':error,'inverse_kick_max_absolute_error':rev,'thermal_error':float(U-thermal)})
    bad=np.r_[q,F,J.ravel(),.001,0.][None,:];bad[0,0]=0;out=np.empty((1,17));assert cpu.mk_batch(1,bad,out)!=0 and gpu.mk_cuda_batch(1,bad,out)!=0
    result={'scope':'Manufactured31positive atoms with non-affine random velocities and prescribed affine acceleration. Component validation, not Kepler or hydro evolution.','rows':rows,'invalid_mass_rejected':True,'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),root/'moment_affine_kick.h',root/'cpu.so',root/'cuda.so')}}
    (root/'atomic_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()

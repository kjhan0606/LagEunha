"""Independent atom oracle and support-criterion controls for material drift."""
import ctypes
import importlib.util
import json
from pathlib import Path
import numpy as np

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
    root=Path(__file__).parent;assert not (root/'atomic_validation.json').exists()
    oracle=load('oracle',Path('/gpfs/kjhan/LagForce/Codex_moment_pressure_kick_20261018/atomic_audit.py'))
    audit=load('audit',root/'audit.py');cpu=ctypes.CDLL(str(root/'cpu.so'));gpu=ctypes.CDLL(str(root/'cuda.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for f in (cpu.md_batch,gpu.md_cuda_batch):f.argtypes=[ctypes.c_size_t,arr,arr];f.restype=ctypes.c_int
    rng=np.random.default_rng(8391);m=rng.uniform(.1,1,41);x=rng.uniform(-.3,.3,(41,2));v=rng.normal(size=x.shape);phi=np.zeros(41)
    q=oracle.moments(m,x,v,3.,phi)[:16];rows=[]
    for dt in (0.,.001,.1,1.,-1.):
        data=np.ascontiguousarray(np.r_[q,dt][None,:]);one=np.empty((1,16));device=np.empty_like(one)
        assert cpu.md_batch(1,data,one)==0 and gpu.md_cuda_batch(1,data,device)==0
        assert np.array_equal(one.view('u8'),device.view('u8'))
        expected=oracle.moments(m,x+dt*v,v,3.,phi)[:16];error=float(np.max(abs(one[0]-expected)));assert error<1e-12
        xx=x+dt*v;lo=xx.min(axis=0)-.1;hi=xx.max(axis=0)+.1;p=np.array([[lo[0],lo[1]],[hi[0],lo[1]],hi,[lo[0],hi[1]]])
        gap,product=audit.support(p,one[0]);assert gap>0 and product>0
        rows.append({'dt':dt,'direct_atom_max_error':error,'support_positive':True})
    p=np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]])
    bad=q.copy();bad[1:3]=q[0]*np.array([3.,0.]);bad[3:6]=q[0]*np.array([9.,0.,0.]);assert audit.support(p,bad)[0]<0
    bad=q.copy();bad[1:3]=0;bad[3:6]=q[0]*np.array([4.,0.,4.]);assert audit.support(p,bad)[1]<0
    (root/'atomic_validation.json').write_text(json.dumps({'rows':rows,'outside_mean_rejected':True,'excess_covariance_rejected':True,'scope':'Manufactured positive atom and impossible-support controls only.'},indent=2)+'\n')

if __name__=='__main__':main()

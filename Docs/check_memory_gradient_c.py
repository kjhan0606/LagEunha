"""Compile and compare C derivative rows to independent SVD on saved geometry.

Read-only geometry test on lageunha. No simulation or live binary replacement.
"""
import ctypes,hashlib,json,subprocess
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

root=Path('/gpfs/kjhan/LagForce/Codex_memory_gradient_probe_20261018')
for label,header in [('old','old_reversible_memory.h'),('new','reversible_memory.h')]:
    source='#include "'+header+'"\nint probe(int n,const double*x,const double*y,int owner,double coverage,double target,double*r,double*w,double*m){return rm_quadratic_row(n,x,y,owner,coverage,target,r,w,m);}\n'
    if label=='new':
        source+='int gradient(int n,const double*x,const double*y,int owner,double coverage,double target,double*r,double*w,double*m,double*gx,double*gy){return rm_quadratic_row_gradient(n,x,y,owner,coverage,target,r,w,m,gx,gy);}\n'
    (root/(label+'.c')).write_text(source)
    subprocess.run(['gcc','-O2','-std=c99','-shared','-fPIC',str(root/(label+'.c')),'-lm','-o',str(root/(label+'.so'))],check=True)
ptr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
signature=[ctypes.c_int,ptr,ptr,ctypes.c_int,ctypes.c_double,ctypes.c_double,ptr,ptr,ptr]
old=ctypes.CDLL(str(root/'old.so'));new=ctypes.CDLL(str(root/'new.so'))
for f in (old.probe,new.probe):f.argtypes=signature;f.restype=ctypes.c_int
new.gradient.argtypes=signature+[ptr,ptr];new.gradient.restype=ctypes.c_int

state=np.loadtxt('/gpfs/kjhan/LagForce/Codex_t6_entropy_force_probes_20261018/adaptive/raw_faces.dat.hgstate.rank0000')
state=state[np.argsort(state[:,0])];x=state[:,3:5]
tree=cKDTree(x%6,boxsize=6)
owners=np.unique(np.r_[np.linspace(0,len(x)-1,256,dtype=int),np.flatnonzero(np.linalg.norm(x-3,axis=1)<.5)[::10]])
maxdiff=0.;maxmoment=0.;maxspin=0.;periodic=0
for owner in owners:
    distance,ids=tree.query(x[owner]%6,k=128)
    dx=x[ids]-x[owner];periodic+=bool(np.any(abs(dx)>3));dx-=6*np.rint(dx/6)
    xx=np.ascontiguousarray(dx[:,0]);yy=np.ascontiguousarray(dx[:,1]);n=len(ids)
    slot=int(np.where(ids==owner)[0][0]);coverage=min(distance[-1],3.)
    rows=[];metas=[]
    for f in (old.probe,new.probe):
        row=np.zeros(n);meta=np.zeros(2);work=np.zeros(7*n)
        assert f(n,xx,yy,slot,coverage,12.,row,work,meta)==0
        rows.append(row);metas.append(meta)
    row=np.zeros(n);meta=np.zeros(2);work=np.zeros(7*n);gx=np.zeros(n);gy=np.zeros(n)
    assert new.gradient(n,xx,yy,slot,coverage,12.,row,work,meta,gx,gy)==0
    assert np.array_equal(rows[0],rows[1]) and np.array_equal(rows[0],row)
    assert np.array_equal(metas[0],metas[1]) and np.array_equal(metas[0],meta)
    radius=meta[0];q=dx/radius
    P=np.column_stack((np.ones(n),q,q[:,0]**2,q[:,0]*q[:,1],q[:,1]**2))
    sqrtw=np.maximum(1-np.sum(q*q,axis=1),0)**2
    u,s,vt=np.linalg.svd(sqrtw[:,None]*P,full_matrices=False)
    for axis,g in enumerate((gx,gy)):
        expected=sqrtw*(u@(vt[:,axis+1]/s))/radius
        maxdiff=max(maxdiff,float(np.max(abs(g-expected))*radius))
        target=np.zeros(6);target[axis+1]=1
        maxmoment=max(maxmoment,float(np.max(abs(P.T@(g*radius)-target))))
    # A manufactured quadratic field with owner vorticity=2*omega.
    omega=2.3
    vx=-omega*yy+.7*xx**2-.2*xx*yy
    vy=omega*xx+.4*yy**2+.3*xx*yy
    maxspin=max(maxspin,abs(.5*(gx@vy-gy@vx)-omega))
assert maxdiff<1e-12 and maxmoment<1e-12 and maxspin<1e-12
result={'rows':len(owners),'periodic_rows':periodic,'old_api_and_gradient_residual_rows_bitwise_equal':True,
        'max_dimensionless_C_vs_SVD_gradient_error':maxdiff,'max_derivative_moment_residual':maxmoment,
        'max_manufactured_spin_error':maxspin,
        'header_sha256':hashlib.sha256((root/'reversible_memory.h').read_bytes()).hexdigest(),
        'scope':'CPU derivative rows on sampled saved t6 geometry. Not CUDA execution, integrated transport, or long-term physics.'}
(root/'analysis.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))

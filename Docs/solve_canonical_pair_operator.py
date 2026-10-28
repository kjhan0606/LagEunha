from pathlib import Path
import json,time,numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsmr
root=Path('/gpfs/kjhan/LagForce/Codex_pair_operator_feasibility_20261017');source=Path('/gpfs/kjhan/LagForce/Codex_gradient_operator_20261017/operator');start=time.monotonic()
a=np.loadtxt(source/'rhs.dat');a=a[np.argsort(a[:,0])];n=len(a);assert n==36864
q=np.fromfile(source/'gradient_operator.bin',dtype='<f8').reshape(-1,6);face=np.loadtxt(source/'raw_faces.dat.rank0000')
pairs=np.concatenate((q[:,:2].astype(np.int64),face[:,:2].astype(np.int64)));pairs=np.sort(pairs,axis=1);pairs=pairs[pairs[:,0]!=pairs[:,1]];keys=np.unique(pairs[:,0]*n+pairs[:,1]);i=keys//n;j=keys%n;m=len(keys)
r=a[j,1:3]-a[i,1:3];r-=6*np.rint(r/6);dx,dy=r.T;d2=dx*dx+dy*dy;assert np.all(d2>0)
fp=np.sort(face[:,:2].astype(np.int64),axis=1);fk=fp[:,0]*n+fp[:,1];idx=np.searchsorted(keys,fk);assert np.all(keys[idx]==fk)
flength=np.hypot(face[:,8]-face[:,6],face[:,9]-face[:,7]);fd=np.hypot(face[:,4],face[:,5]);s0=np.zeros(m);counts=np.zeros(m);np.add.at(s0,idx,flength/fd);np.add.at(counts,idx,1);s0=np.divide(s0,counts,out=np.zeros_like(s0),where=counts>0)
vol=a[:,6];h=np.sqrt(vol);edge=np.arange(m);rr=[];cc=[];vv=[]
for node,sign in [(i,1),(j,-1)]:
 terms=[sign*dx/h[node],sign*dy/h[node],dx*dx/vol[node],dx*dy/vol[node],dy*dy/vol[node]]
 for k,values in enumerate(terms):rr.append(5*node+k);cc.append(edge);vv.append(values)
A=coo_matrix((np.concatenate(vv),(np.concatenate(rr),np.concatenate(cc))),shape=(5*n,m)).tocsr();del rr,cc,vv,q,pairs
b=np.zeros((n,5));b[:,2]=b[:,4]=2;b=b.ravel();initial=A@s0-b
print(json.dumps({'particles':n,'edges':m,'nonzeros':A.nnz,'initial_max_residual':float(abs(initial).max())}),flush=True)
result=lsmr(A,-initial,atol=1e-11,btol=1e-11,conlim=1e12,maxiter=1000,show=True);s=s0+result[0];residual=(A@s-b).reshape(n,5)
np.savez(root/'coefficients.npz',i=i,j=j,dx=dx,dy=dy,s=s,s0=s0,volume=vol)
dv=a[j,3:5]-a[i,3:5];pair=.5*s*np.sum(r*dv,axis=1);vdot=np.bincount(i,weights=pair,minlength=n)+np.bincount(j,weights=pair,minlength=n)
closure=np.column_stack((np.bincount(i,weights=s*dx,minlength=n)-np.bincount(j,weights=s*dx,minlength=n),np.bincount(i,weights=s*dy,minlength=n)-np.bincount(j,weights=s*dy,minlength=n)))
radius=np.hypot(a[:,1]-3,a[:,2]-3);regions={};p0=2.7777796321447566e-5
for name,sel in [('inner',radius<.5),('disk_inner_edge',(radius>=.5)&(radius<.6)),('disk_core',(radius>=.6)&(radius<1.9))]:
 acc=p0*np.linalg.norm(closure[sel],axis=1)/a[sel,7];weights=a[sel,7]
 regions[name]={'max_closure_scaled_residual':float(np.max(np.abs(residual[sel,:2]))),'max_affine_moment_residual':float(np.max(np.abs(residual[sel,2:]))),'constant_pressure_acceleration_RMS':float(np.sqrt(np.average(acc**2,weights=weights))),'circular_divergence_RMS':float(np.sqrt(np.average((vdot[sel]/vol[sel])**2,weights=weights))),'original_geometric_divergence_RMS':float(np.sqrt(np.average((a[sel,12]/vol[sel])**2,weights=weights)))}
out={'scope':'Static feasibility on actual full canonical t2 graph, union of C QR supports and Voronoi faces. B_i(v)=0.5 sum_j s_ij r_ij dot(v_j-v_i), symmetric scalar s. Constraints sum_j s r=0 and sum_j s r r^T=2 V_i I. Not a time-evolution solver.','n':n,'edges':m,'lsmr_stop_code':int(result[1]),'iterations':int(result[2]),'norm_residual':float(result[3]),'max_constraint_residual':float(np.max(np.abs(residual))),'coefficient_min':float(s.min()),'coefficient_max':float(s.max()),'negative_fraction':float(np.mean(s<0)),'regions':regions,'seconds':time.monotonic()-start,'converged_to_requested_constraint_tolerance':bool(np.max(np.abs(residual))<1e-8)}
(root/'result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)

from pathlib import Path
import numpy as np,json
root=Path(__file__).parent;source=Path('/gpfs/kjhan/LagForce/Codex_gradient_operator_20261017/operator')
z=np.load(root/'response.npz');i,j=z['i'],z['j'];n=36864
f=np.loadtxt(source/'raw_faces.dat.common.rank0000');keys=f[:,0].astype(int)*n+f[:,1].astype(int);order=np.argsort(keys);f=f[order];keys=keys[order];idx=np.searchsorted(keys,i*n+j);assert np.array_equal(keys[idx],i*n+j)
f=f[idx];A=f[:,4];d=np.linalg.norm(f[:,2:4],axis=1);required=np.divide(2*d*z['required_face_flux_correction'],A,out=np.zeros_like(A),where=A>0)
a=np.loadtxt(source/'rhs.dat');a=a[np.argsort(a[:,0])];r=np.linalg.norm(a[:,1:3]-3,axis=1)
lookup={(int(ii),int(jj)):k for k,(ii,jj) in enumerate(zip(i,j)) if A[k]>0}
neighbors=[set() for _ in range(n)]
for ii,jj in lookup:neighbors[ii].add(jj);neighbors[jj].add(ii)
results={}
for name,mask in [('inner',r<.5),('disk_core',(r>=.6)&(r<1.9))]:
 cycles=[]
 for ii in np.flatnonzero(mask):
  for jj in neighbors[ii]:
   if jj<=ii or not mask[jj]:continue
   for kk in neighbors[ii]&neighbors[jj]:
    if kk<=jj or not mask[kk]:continue
    ij,jk,ik=lookup[(ii,jj)],lookup[(jj,kk)],lookup[(ii,kk)]
    terms=[required[ij],required[jk],-required[ik]];res=sum(terms)
    cycles.append((abs(res),int(ii),int(jj),int(kk),res,terms,float(np.mean(d[[ij,jk,ik]]**2))))
 worst=max(cycles);results[name]={'triangle_count':len(cycles),'cycle_sum_RMS':float(np.sqrt(np.mean([x[4]**2 for x in cycles]))),'max_abs_cycle_sum':worst[0],'example_ids':list(worst[1:4]),'example_signed_required_weight_differences':worst[5],'example_signed_cycle_sum':worst[4],'example_cycle_sum_over_mean_edge_d_squared':worst[4]/worst[6]}
 assert worst[0]>1e-10
out={'scope':'Interior 3-edge cycle certificate, independent of frozen-weight boundary choice. A scalar weight potential has identically zero sum of differences around every cycle. Nonzero required sum rules out exact instantaneous matching of all reconstructed face normal velocities at this geometry. Does not rule out matching volumes or another reconstruction.','regions':results}
(root/'cycle_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

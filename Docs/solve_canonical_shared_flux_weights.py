from pathlib import Path
import numpy as np,json,hashlib,re
from scipy.sparse import coo_matrix
root=Path(__file__).parent;source=Path('/gpfs/kjhan/LagForce/Codex_gradient_operator_20261017/operator')
a=np.loadtxt(source/'rhs.dat');a=a[np.argsort(a[:,0])];n=len(a)
q=np.fromfile(source/'gradient_operator.bin',dtype='<f8').reshape(-1,6);qi,qj=q[:,:2].astype(int).T
Wx=coo_matrix((q[:,2],(qi,qj)),shape=(n,n)).tocsr();Wy=coo_matrix((q[:,3],(qi,qj)),shape=(n,n)).tocsr()
f=np.loadtxt(source/'raw_faces.dat.common.rank0000');fi,fj=f[:,:2].astype(int).T
# Use one i<j face; verify every retained edge has its reciprocal before using it.
keys=fi*n+fj;order=np.argsort(keys);keys=keys[order];f=f[order];fi,fj=f[:,:2].astype(int).T
rev=np.searchsorted(keys,fj*n+fi);paired=(rev<len(keys));paired[paired]=keys[rev[paired]]==fj[paired]*n+fi[paired]
unpaired=f[~paired,:2].astype(int).tolist();sel=(fi<fj)&paired;g=f[sel];back=f[rev[sel]];i,j=g[:,:2].astype(int).T
r=g[:,2:4];d=np.linalg.norm(r,axis=1);e=r/d[:,None];A=g[:,4];ri=g[:,5:7];rj=ri-r
backA=back[:,4]
assert np.max(np.abs(A-backA))<1e-10
assert np.max(np.abs(r+back[:,2:4]))<1e-10
assert np.max(np.abs(rj-back[:,5:7]))<1e-10
Ae=A[:,None]*e

def scatter(values,idx):
 z=np.zeros((n,)+values.shape[1:]);np.add.at(z,idx,values);return z

def B(v):
 gx=Wx@v;gy=Wy@v
 vf=.5*(v[i]+v[j]+gx[i]*ri[:,0,None]+gy[i]*ri[:,1,None]+gx[j]*rj[:,0,None]+gy[j]*rj[:,1,None])
 flux=np.sum(Ae*vf,axis=1)
 return scatter(flux,i)-scatter(flux,j)

def BT(p):
 z=.5*(p[i]-p[j])[:,None]*Ae
 direct=scatter(z,i)+scatter(z,j)
 sx=scatter(z*ri[:,0,None],i)+scatter(z*rj[:,0,None],j)
 sy=scatter(z*ri[:,1,None],i)+scatter(z*rj[:,1,None],j)
 return direct+Wx.T@sx+Wy.T@sy


from scipy.sparse.linalg import spsolve
import time
started=time.monotonic()
xy=a[:,1:3]-3;radii=np.linalg.norm(xy,axis=1);active=radii<2.5
v=a[:,3:5];gx=Wx@v;gy=Wy@v
vf=.5*(v[i]+v[j]+gx[i]*ri[:,0,None]+gy[i]*ri[:,1,None]+gx[j]*rj[:,0,None]+gy[j]*rj[:,1,None])
flux=np.sum(Ae*vf,axis=1)
t=ri-.5*r;t-=np.sum(t*e,axis=1)[:,None]*e
geomflux=A*(.5*np.sum((v[i]+v[j])*e,axis=1)+np.sum((v[i]-v[j])*t,axis=1)/d)
delta=flux-geomflux
conductance=A/(2*d)
L=coo_matrix((np.concatenate((conductance,conductance,-conductance,-conductance)),(np.concatenate((i,j,i,j)),np.concatenate((i,j,j,i)))),shape=(n,n)).tocsr()
rhs=scatter(delta,i)-scatter(delta,j)
ids=np.flatnonzero(active);La=L[ids][:,ids].tocsc();wd=np.zeros(n);wd[ids]=spsolve(La,rhs[ids])
assert np.isfinite(wd).all()
volume_residual=L@wd-rhs
relative=np.linalg.norm(volume_residual[active])/np.linalg.norm(rhs[active]);assert relative<1e-9
weightflux=conductance*(wd[i]-wd[j]);residual_flux=delta-weightflux
# Weighted least-squares face-speed fit: sum A*d*(delta_normal-delta_w/(2d))^2.
# Euler equation is the same volume Laplacian; residual has zero divergence in active cells.
velres=np.divide(residual_flux,A,out=np.zeros_like(A),where=A>0)
veltarget=np.divide(delta,A,out=np.zeros_like(A),where=A>0)
regions={}
for name,mask in [('inner',radii<.5),('disk_edge',(radii>=.5)&(radii<.6)),('disk_core',(radii>=.6)&(radii<1.9)),('active',active)]:
 edge=mask[i]&mask[j]&(A>0);w=A[edge]*d[edge]
 wrms=lambda z:float(np.sqrt(np.average(z[edge]**2,weights=w)))
 regions[name]={'cells':int(mask.sum()),'internal_edges':int(edge.sum()),'weight_rate_over_volume_RMS':float(np.sqrt(np.average((wd[mask]/a[mask,6])**2,weights=a[mask,7]))),'weight_rate_over_volume_max':float(np.max(abs(wd[mask]/a[mask,6]))),'required_face_normal_correction_RMS':wrms(veltarget),'unrealized_face_normal_correction_RMS':wrms(velres),'unrealized_to_required_RMS_ratio':wrms(velres)/wrms(veltarget),'volume_rate_fit_relative_L1':float(np.sum(abs(volume_residual[mask]))/np.sum(abs(rhs[mask])))}
# Linearized edge-weight metric at selected dt; not a positivity or future-growth proof.
edge=active[i]|active[j]
diff_rate=(wd[i]-wd[j])/d**2
result={'source':str(source),'scope':'Instantaneous W=0 linear response at t2 circular geometry; fixed generator velocity and mass. Unknown W is power-diagram squared-length weight, not its square root. Frozen r>=2.5 weights held fixed. Shared paired faces only; cap contributions outside active region omitted.','equation':'L_ii=sum A/(2d), L_ij=-A/(2d); L Wdot=B_shared v-B_geom v on active cells.','active_linear_solve_relative_L2_residual':float(relative),'regions':regions,'weight_difference_rate_over_d_squared_max_touching_active':float(abs(diff_rate[edge]).max()),'inferred_max_edge_weight_difference_over_d_squared_at_dt0025':float(.0025*abs(diff_rate[edge]).max()),'frozen_induced_total_volume_rate':float((L@wd)[~active].sum()),'active_target_total_volume_correction':float(rhs[active].sum()),'energy_work_required_if_old_pressure_force_retained':float(-np.dot(a[active,5],rhs[active])),'seconds':time.monotonic()-started,'limitations':['Only instantaneous linearized volume feasibility; does not establish nonlinear Laguerre-cell positivity or bounded weights over time.','Residual face-speed mismatch can circulate even when volume constraints are exactly satisfied. Matching cell volumes does not prove material face motion.','Old pressure force plus new volume rate breaks energy work unless the matched transpose force or another explicitly consistent energy model is used.','No production source or evolving run changed.']}
np.savez(root/'response.npz',weight_rate=wd,required_face_flux_correction=delta,weight_face_flux_correction=weightflux,residual_face_flux=residual_flux,volume_rhs=rhs,i=i,j=j)
(root/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

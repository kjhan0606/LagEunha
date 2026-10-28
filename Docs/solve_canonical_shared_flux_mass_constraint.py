from pathlib import Path
import numpy as np,json,time
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
root=Path(__file__).parent;source=Path('/gpfs/kjhan/LagForce/Codex_gradient_operator_20261017/operator');prior=Path('/gpfs/kjhan/LagForce/Codex_shared_flux_weight_response_20261017')
a=np.loadtxt(source/'rhs.dat');a=a[np.argsort(a[:,0])];n=len(a);z=np.load(prior/'response.npz');i,j=z['i'],z['j'];delta=z['required_face_flux_correction'];wv=z['weight_rate']
f=np.loadtxt(source/'raw_faces.dat.common.rank0000');keys=f[:,0].astype(int)*n+f[:,1].astype(int);order=np.argsort(keys);keys=keys[order];f=f[order];idx=np.searchsorted(keys,i*n+j);assert np.array_equal(keys[idx],i*n+j);f=f[idx]
A=f[:,4];dist=np.linalg.norm(f[:,2:4],axis=1);c=A/(2*dist);assert np.all(c>0)
r=np.linalg.norm(a[:,1:3]-3,axis=1);active=r<2.5;ids=np.flatnonzero(active);rho=a[:,7]/a[:,6]
def div(x):
 if x.ndim==1:return np.bincount(i,weights=x,minlength=n)-np.bincount(j,weights=x,minlength=n)
 out=np.zeros((n,x.shape[1]));np.add.at(out,i,x);np.add.at(out,j,-x);return out
def evaluate(w):
 q=delta-c*(w[i]-w[j]);up=np.where(q>=0,rho[i],rho[j]);mf=up*q
 return q,up,mf,div(mf),float(np.sum(.5*up*q*q/c))
def matrix(up):
 t=c*up
 return coo_matrix((np.concatenate((t,t,-t,-t)),(np.concatenate((i,j,i,j)),np.concatenate((i,j,j,i)))),shape=(n,n)).tocsr()[ids][:,ids].tocsc()
started=time.monotonic();w=wv.copy();scale=np.linalg.norm(evaluate(np.zeros(n))[3][active]);history=[]
for iteration in range(50):
 q,up,mf,res,obj=evaluate(w);rel=np.linalg.norm(res[active])/scale
 history.append({'iteration':iteration,'relative_mass_balance_L2':float(rel),'objective':obj})
 if rel<1e-11:break
 dw=np.zeros(n);dw[ids]=spsolve(matrix(up),res[active]);slope=-np.dot(res[active],dw[active])
 assert slope<0
 alpha=1.
 while alpha>2**-24:
  trial=w+alpha*dw;value=evaluate(trial)[4]
  if value<=obj+1e-4*alpha*slope:break
  alpha*=.5
 assert alpha>2**-24,'line search failed'
 history[-1]['alpha']=alpha;w=trial
else:raise RuntimeError('mass solve did not converge')
assert rel<1e-11
qv,upv,mfv,resv,_=evaluate(wv)
volume_delta=div(c*(w[i]-w[j]))-div(c*(wv[i]-wv[j]))
regions={}
for name,sel in [('inner',r<.5),('disk_edge',(r>=.5)&(r<.6)),('disk_core',(r>=.6)&(r<1.9)),('active',active)]:
 rms=lambda v:float(np.sqrt(np.average(v[sel]**2,weights=a[sel,7])))
 edge=sel[i]&sel[j];dv=a[i,3:5]-a[j,3:5]
 regions[name]={'n':int(sel.sum()),'volume_matching_fractional_mass_rate_RMS':rms(resv/a[:,7]),'volume_matching_fractional_mass_rate_max':float(np.max(abs(resv[sel]/a[sel,7]))),'mass_matching_fractional_mass_rate_RMS':rms(res/a[:,7]),'mass_matching_fractional_mass_rate_max':float(np.max(abs(res[sel]/a[sel,7]))),'additional_volume_divergence_RMS':rms(volume_delta/a[:,6]),'weight_rate_over_volume_RMS':rms(w/a[:,6]),'internal_face_mass_circulation_L1':float(np.sum(abs(mf[edge]))),'first_order_upwind_velocity_mixing_power':float(.5*np.sum(abs(mf[edge])*np.sum(dv[edge]**2,axis=1)))}
# Donor transport alone: exact paired mass/momentum/total-energy cancellation.
vel=a[:,3:5];donor=np.where(q>=0,i,j);mdot=-div(mf);pdot=-div(mf[:,None]*vel[donor]);specific_energy=.5*np.sum(vel*vel,axis=1)+a[:,5]*a[:,6]/(a[:,7]*(1.6666669845581055-1));edot=-div(mf*specific_energy[donor])
kinetic_rate=float(np.sum(vel*pdot)-.5*np.dot(np.sum(vel*vel,axis=1),mdot));mixing=-float(.5*np.sum(abs(mf)*np.sum((vel[i]-vel[j])**2,axis=1)))
assert abs(kinetic_rate-mixing)<1e-12*max(abs(mixing),1e-30)
out={'scope':'Static diagnostic of scalar-weight rates enforcing zero NET mass transfer per active cell, with piecewise constant upwind density on shared reconstructed face-velocity residuals. Not a production ALE/Riemann solver or long-time result.','equation':'sum_j rho_upwind*[delta_face_volume_flux-A/(2d)*(Wdot_i-Wdot_j)] = 0, frozen weight rates=0.','iterations':history,'seconds':time.monotonic()-started,'regions':regions,'transport_algebra':{'global_mass_rate':float(mdot.sum()),'global_momentum_rate':pdot.sum(axis=0).tolist(),'global_total_energy_rate':float(edot.sum()),'kinetic_transport_rate':kinetic_rate,'negative_upwind_mixing_power':mixing,'active_mass_rate_sum':float(mdot[active].sum()),'frozen_mass_rate_L1':float(abs(mdot[~active]).sum())},'limitations':['Net zero cell mass does not imply material faces or zero exchange through individual faces.','Remaining momentum and total-energy fluxes must be retained; setting them to zero is not an ALE discretization.','First-order donor mixing is only a diagnostic, not a recommended shear-stable transport method.','No nonlinear mesh update, time integration, shock test, positivity test, boundary model, or long-term weight bound is established.']}
np.savez(root/'response.npz',weight_rate=w,face_mass_flux=mf,face_relative_volume_flux=q,active=active,i=i,j=j)
(root/'result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

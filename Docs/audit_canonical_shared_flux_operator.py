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

xy=a[:,1:3]-3;radii=np.linalg.norm(xy,axis=1);active=radii<2.5
rng=np.random.default_rng(20261017);vel=rng.normal(size=(n,2));pressure=rng.uniform(.5,1.5,n)
force=BT(pressure);volume=B(vel);scale=np.sum(abs(vel*force))+np.sum(abs(pressure*volume))
energy_res=abs(np.sum(vel*force)-np.dot(pressure,volume))/scale
assert energy_res<1e-12
pactive=np.where(active,pressure,0);fa=BT(pactive)
rotation=np.column_stack((-xy[:,1],xy[:,0]));tx=np.tile([1.,0.],(n,1));ty=np.tile([0.,1.],(n,1))
checks={'translation_x_max_active':float(abs(B(tx)[active]).max()),'translation_y_max_active':float(abs(B(ty)[active]).max()),'rotation_max_active':float(abs(B(rotation)[active]).max()),'isotropic_expansion_relative_max_active':float(abs((B(xy)-2*a[:,6])[active]/a[active,6]).max()),'global_volume_sum_random_velocity':float(volume.sum()),'constant_pressure_force_max':float(abs(BT(np.ones(n))).max()),'adjoint_energy_relative_residual':float(energy_res),'net_force_active_pressure':fa.sum(axis=0).tolist(),'net_torque_active_pressure':float(np.sum(xy[:,0]*fa[:,1]-xy[:,1]*fa[:,0]))}
affine={}
for aa in range(2):
 for bb in range(2):
  vtest=np.zeros((n,2));vtest[:,aa]=xy[:,bb];expected=a[:,6] if aa==bb else np.zeros(n)
  affine[f'dv{aa}_dx{bb}']=float(np.max(np.abs((B(vtest)-expected)[active]/a[active,6])))
checks['all_affine_relative_max_active']=affine
assert max(affine.values())<1e-10
assert checks['constant_pressure_force_max']==0
assert checks['rotation_max_active']<1e-12
vd=B(a[:,3:5]);regions={} 
for name,s in [('inner',radii<.5),('disk_edge',(radii>=.5)&(radii<.6)),('disk_core',(radii>=.6)&(radii<1.9))]:
 rms=lambda z:float(np.sqrt(np.average(z[s]**2,weights=a[s,7])))
 regions[name]={'n':int(s.sum()),'candidate_divergence_rms':rms(vd/a[:,6]),'original_geometric_divergence_rms':rms(a[:,12]/a[:,6]),'raw_QR_divergence_rms':rms(Wx@a[:,3]+Wy@a[:,4])}
gamma=float(np.float32(float(re.search(r'define Gamma_hydro\s*=\s*([\d.eE+-]+)',(source/'params.dat').read_text())[1])));p0=2.7777796321447566e-5
for name,mask in [('inner',radii<.5),('disk_edge',(radii>=.5)&(radii<.6)),('disk_core',(radii>=.6)&(radii<1.9))]:
 w=a[mask,7];pn=a[mask,5]/p0;fluct=pn-np.average(pn,weights=w)
 for label,rate in [('candidate',vd),('original',a[:,12])]:
  dp=-gamma*pn*rate[mask]/a[mask,6]
  regions[name][label+'_compression_pressure_variance_rate']=float(2*np.average(fluct*dp,weights=w))
result={'scope':'Static feasibility of shared reconstructed face-volume flux using actual C QR gradients and paired real faces at t2 circular diagnostic geometry. No time evolution, no production physics changed.','source':str(source),'operator_sha256':hashlib.sha256((source/'gradient_operator.bin').read_bytes()).hexdigest(),'gamma':gamma,'edges':len(i),'unpaired_excluded':unpaired,'checks':checks,'regions':regions,'limitations':['Far frozen cap edges are absent; this operator is not a complete boundary model.','Stored geometric volumes remain the normalization here; adopting independent volumes requires positivity and physical meaning.','Weight-based realization would need a compatible work/force model and proof against secular weight growth.','Affine and conservation checks do not establish shock or long-term Kepler stability.']}
(root/'audit_final.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

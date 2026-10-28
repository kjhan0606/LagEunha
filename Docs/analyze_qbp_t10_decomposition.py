"""Actual QR map at fixed evolved geometry; circular counterfactual velocity."""
import json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from audit_canonical_pair_history import gravity_coefficient
p=Path('/gpfs/kjhan/LagForce/Codex_qbp_t10_decomposition_20261018')
assert (p/'exit.code').read_text().strip()=='0'
s=np.loadtxt(p/'raw_faces.dat.rhs.rank0000');s=s[np.argsort(s[:,0])]
c=np.loadtxt(p/'raw_faces.dat.cache.rank0000');c=c[np.argsort(c[:,0])]
assert np.all(c[:,5]==2) and np.all(c[:,6]==0),'requires valid non-arc fits'
f=np.loadtxt(p/'raw_faces.dat.rhs.faces.rank0000')
op=np.fromfile(p/'gradient_operator.bin').reshape(-1,6);i,j=op[:,:2].astype(int).T;n=len(s)
def gradient(v):
 return np.column_stack([np.bincount(i,weights=op[:,k]*v[j,a],minlength=n) for a,k in [(0,2),(0,3),(1,2),(1,3)]])
x=s[:,1:3]-3;r=np.linalg.norm(x,axis=1);v=s[:,3:5]
vc=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack((-x[:,1],x[:,0]))
ii,jj=f[:,:2].astype(int).T;e=f[:,6:8];d=f[:,5]
def signals(v):
 g=gradient(v);du=np.sum((v[ii]-v[jj])*e,axis=1);gs=g[ii]+g[jj]
 delta=du+.5*d*(gs[:,0]*e[:,0]**2+(gs[:,1]+gs[:,2])*e[:,0]*e[:,1]+gs[:,3]*e[:,1]**2)
 sf=np.sqrt(g[:,0]**2+g[:,3]**2+.5*(g[:,1]+g[:,2])**2)
 chi=np.minimum(1,np.maximum(-(g[:,0]+g[:,3]),0)/np.maximum(sf,1e-300))
 cf=np.maximum(chi[ii],chi[jj]);qin=(1-cf)*delta+cf*du
 return g,du,delta,cf,qin
actual=signals(v);circ=signals(vc);pert=signals(v-vc)
def rel(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-300))
out={'scope':'Frozen t10 positions and actual exported QR coefficients. Circular and perturbation delta add linearly; recomputed chi/minmod/heat do not. No evolved counterfactual trajectory.',
 'gradient_reproduction_max_abs':float(np.max(abs(actual[0]-c[:,7:11]))),
 'recorded_delta_relative_l2':rel(actual[2],f[:,19]),'recorded_qin_relative_l2':rel(actual[4],f[:,18]),
 'linear_delta_decomposition_max_abs':float(np.max(abs(actual[2]-circ[2]-pert[2]))),'regions':{}}
assert out['gradient_reproduction_max_abs']<1e-10
assert out['recorded_qin_relative_l2']<1e-4
for name,mask in [('inner',r[ii]<.5),('edge',(r[ii]>=.5)&(r[ii]<.6)),('core',(r[ii]>=.6)&(r[ii]<1.9))]:
 w=f[mask,4];rms=lambda a:float(np.sqrt(np.average(a[mask]**2,weights=w)))
 out['regions'][name]={'directed_faces':int(mask.sum()),'area_weighted_delta_RMS':{'actual':rms(actual[2]),'circular':rms(circ[2]),'perturbation':rms(pert[2])},
 'area_weighted_qin_RMS':{'actual':rms(actual[4]),'circular_recomputed_chi':rms(circ[4])},
 'area_weighted_chi_mean':{'actual':float(np.average(actual[3][mask],weights=w)),'circular':float(np.average(circ[3][mask],weights=w))}}
with (p/'analysis.json').open('x') as outfp:json.dump(out,outfp,indent=2)
print(json.dumps(out,indent=2))

"""TCF pair identities, using recorded states and synthetic tangential offsets.

No face centroids were recorded in this dump: these offsets are deliberately
synthetic. This tests algebra, not a replay of actual cell force assembly.
"""
import json
from pathlib import Path
import numpy as np

p=Path('/gpfs/kjhan/LagForce/Codex_TCF_rhs_identity_20261017/t6')
f=np.loadtxt(p/'rhs.dat.faces'); s=np.loadtxt(p/'rhs.dat')
s=s[np.argsort(s[:,0])]
ids=s[:,0].astype(int)
ii=np.searchsorted(ids,f[:,0].astype(int)); jj=np.searchsorted(ids,f[:,1].astype(int))
assert np.array_equal(ids[ii],f[:,0]) and np.array_equal(ids[jj],f[:,1])
A=f[:,4]; d=f[:,5]; e=f[:,6:8]
vi=s[ii,3:5]; vj=s[jj,3:5]; pi=s[ii,5]; pj=s[jj,5]; Q=f[:,14]
# Re-normalize rounded normals so constructed offsets are truly tangential.
e=e/np.linalg.norm(e,axis=1)[:,None]
t=.37*d[:,None]*np.column_stack((-e[:,1],e[:,0]))
def dot(x,y): return np.sum(x*y,axis=1)
def assemble(vi,vj,pi,pj,e,t):
    ui=dot(vi,e); uj=dot(vj,e)
    F=-A[:,None]*((.5*(pi+pj)+Q)[:,None]*e+((pj-pi)/d)[:,None]*t)
    E=-A*(.5*(pj*ui+pi*uj)+(pj*dot(vi,t)-pi*dot(vj,t))/d+.5*Q*(ui+uj))
    V=A*(.5*(ui+uj)+dot(vi-vj,t)/d)
    H=.5*A*Q*(ui-uj)
    return F,E,V,H
F,E,V,H=assemble(vi,vj,pi,pj,e,t)
Fr,Er,Vr,Hr=assemble(vj,vi,pj,pi,-e,t)
b=np.array([7.,-3.]); Fb,Eb,Vb,Hb=assemble(vi+b,vj+b,pi,pj,e,t)
def maximum(x): return float(np.max(abs(x)))
out={
 'rows':len(f),'tangent_scope':'synthetic t=0.37*d*perpendicular(e); not recorded geometry',
 'force_reversal_max_abs':maximum(F+Fr),
 'energy_reversal_max_abs':maximum(E+Er),
 'volume_reversal_max_abs':maximum(V+Vr),
 'heat_reversal_max_abs':maximum(H-Hr),
 'boost_force_max_abs':maximum(Fb-F),
 'boost_energy_covariance_max_abs':maximum(Eb-E-dot(F,np.broadcast_to(b,F.shape))),
 'boost_volume_covariance_max_abs':maximum(Vb-V-A*dot(e,np.broadcast_to(b,e.shape))),
 'boost_heat_max_abs':maximum(Hb-H),
 'face_thermal_identity_max_abs':maximum(E-dot(vi,F)+pi*V-H-A*pi*dot(vi,e)),
 'identity':'E_i-v_i.F_i+P_i Vdot_i-H_i = A P_i v_i.e; last term cancels only on a closed cell',
 'production_dump_closure_thermal_defect_max_abs':maximum(s[:,5]*(s[:,3]*s[:,15]+s[:,4]*s[:,16])),
 'limitations':'No reconstructed geometry, QBP, MPI, CUDA, HG, or RK integration tested here.'
}
radius=np.hypot(s[:,1]-3,s[:,2]-3)
defect=s[:,5]*(s[:,3]*s[:,15]+s[:,4]*s[:,16])
out['coordinate_center']=[3.,3.]
out['closure_defect_by_radius']={name:{'count':int(mask.sum()),
    'max_abs':maximum(defect[mask]),'rms':float(np.sqrt(np.mean(defect[mask]**2))),
    'sumAn_max_norm':float(np.max(np.linalg.norm(s[mask,13:15],axis=1))),
    'sumAe_max_norm':float(np.max(np.linalg.norm(s[mask,15:17],axis=1)))}
    for name,mask in [('r_lt_1p9',radius<1.9),('r_ge_1p9',radius>=1.9)]}
out['largest_closure_cells']=[{'id':int(s[k,0]),'radius':float(radius[k]),
    'thermal_closure_term':float(defect[k])} for k in np.argsort(abs(defect))[-10:][::-1]]
root=Path('/gpfs/kjhan/LagForce/Codex_tcf_pair_identities_corrected_20261018')
root.mkdir(exist_ok=False)
(root/'result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))

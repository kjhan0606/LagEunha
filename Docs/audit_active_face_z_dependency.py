"""Read saved TCF faces: distinguish overwritten legacy impedance from applied FACE_Z."""
from pathlib import Path
import json,numpy as np
base=Path('/gpfs/kjhan/LagForce');paths={'t2':base/'Codex_t2_onset_rhs_20261017/t2','t4':base/'Codex_TCF_rhs_identity_20261017/t4','t6':base/'Codex_TCF_rhs_identity_20261017/t6','t10':base/'Codex_TCF_rhs_t10_20261017/t10'};rows=[]
for name,p in paths.items():
 f=np.loadtxt(p/'rhs.dat.faces');assert f.shape[1]==21 and np.isfinite(f).all()
 A,du,Zold,Q,Zused,qin=(f[:,j] for j in (4,9,10,14,15,18));mm=np.where(qin*du>0,np.sign(qin)*np.minimum(abs(qin),abs(du)),0.);pred=Zused*mm;assert np.allclose(Q,pred,rtol=3e-5,atol=1e-16)
 heat=.5*A*Q*du;legacy_heat=.5*A*Zold*mm*du;active=abs(Q)>1e-20;different=abs(Zold-Zused)>1e-5*np.maximum(np.maximum(abs(Zold),abs(Zused)),1e-30)
 rows.append({'case':name,'directed_faces':len(f),'nonzero_Q_faces':int(active.sum()),'Z_difference_above_relative_1e_minus5':int(different.sum()),'active_Q_faces_with_Z_difference':int((active&different).sum()),'Q_identity_relative_L2_error':float(np.linalg.norm(pred-Q)/max(np.linalg.norm(Q),1e-300)),'total_directed_heat':float(heat.sum()),'substituting_legacy_Z_heat_ratio':float(legacy_heat.sum()/heat.sum()),'min_directed_heat':float(heat.min())})
out={'scope':__doc__,'rows':rows,'limits':['Printed faces ~7significant digits.','Legacy-Z substitution is a frozen diagnostic, not an evolved alternative or recommended correction.','Pre-HG saved hydro face work, not complete step energy accounting.']}
root=base/'Codex_active_tcf_pressure_20261018'
with (root/'face_z_dependency.json').open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
print(json.dumps(out,indent=2))

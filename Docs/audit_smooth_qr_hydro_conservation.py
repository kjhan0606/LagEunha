"""Raw TCF total-energy and thermodynamic identities at preserved states."""
import json,argparse
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root',type=Path,default=Path('/gpfs/kjhan/LagForce/Codex_smooth_qr_probes_20261018'))
root=parser.parse_args().root;out={}
for name in ['off_plus','on_plus','on_minus']:
 d=root/name;a=np.loadtxt(d/'raw_faces.dat.rhs.rank0000');f=np.loadtxt(d/'raw_faces.dat.rhs.faces.rank0000')
 force=a[:,8:10];power=np.sum(a[:,3:5]*force,axis=1)
 # Column10 is dte, the TOTAL energy rate, not the internal energy rate.
 total_rate=a[:,10];internal_rate=total_rate-power
 reversible=-a[:,5]*a[:,12];heating=a[:,11]
 residual=internal_rate-reversible-heating
 mask=np.hypot(a[:,1]-3,a[:,2]-3)<1.9
 scale=abs(internal_rate)+abs(reversible)+abs(heating)
 rel=float(abs(residual[mask]).sum()/scale[mask].sum())
 momentum=float(np.linalg.norm(force.sum(axis=0))/np.linalg.norm(force,axis=1).sum())
 energy=float(abs(total_rate.sum())/abs(total_rate).sum())
 assert rel<1e-9 and momentum<1e-10 and energy<1e-10
 entropy=f[:,14]*f[:,9];assert np.all(entropy>=0)
 out[name]={'force_sum_over_sum_magnitudes':momentum,'total_energy_rate_sum':float(total_rate.sum()),
  'total_energy_rate_relative_residual':energy,'interior_thermodynamic_identity_relative_L1_residual':rel,
  'face_count':len(f),'negative_Q_du_faces':int(np.count_nonzero(entropy<0))}
result={'cases':out,'supersedes':'hydro_conservation.json: exploratory check mislabeled total-energy dte as internal energy and double-counted kinetic work. Preserve that file as invalid, do not use its energy residual.',
 'scope':'Raw TCF RHS before gravity/HG/frozen-cell evolution constraints. Checks sum(dE)=0 and dE-v.F=-P*dV+TdS; face heating sign from rounded records. Not an accepted-step or long-time conservation claim.'}
with (root/'hydro_conservation_corrected.json').open('x') as f:json.dump(result,f,indent=2)
print(json.dumps(result,indent=2))

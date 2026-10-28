"""Check saved TCF closure diagnostics with the verified (3,3) center."""
import json
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
cases={'t2':'Codex_t2_onset_rhs_20261017/t2',
       't4':'Codex_TCF_rhs_identity_20261017/t4',
       't6':'Codex_TCF_rhs_identity_20261017/t6',
       't10':'Codex_TCF_rhs_t10_20261017/t10'}
result={'center':[3,3],'scope':'Saved case labels; snapshot times not independently re-read.', 'cases':{}}
for label,path in cases.items():
    s=np.loadtxt(base/path/'rhs.dat')
    r=np.hypot(s[:,1]-3,s[:,2]-3)
    z=s[:,5]*np.sum(s[:,3:5]*s[:,15:17],axis=1)
    result['cases'][label]={}
    for name,mask in [('inner',r<.5),('disk', (r>=.5)&(r<1.9)),('outer',r>=1.9)]:
        t=s[mask]; d=z[mask]
        result['cases'][label][name]={
            'cells':len(t),'closure_term_max_abs':float(np.max(abs(d))),
            'closure_term_rms':float(np.sqrt(np.mean(d*d))),
            'heat_rms':float(np.sqrt(np.mean(t[:,11]**2))),
            'sumAn_max_norm':float(np.max(np.linalg.norm(t[:,13:15],axis=1))),
            'sumAe_max_norm':float(np.max(np.linalg.norm(t[:,15:17],axis=1)))}
p=base/'Codex_tcf_closure_history_20261018'; p.mkdir(exist_ok=False)
(p/'result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

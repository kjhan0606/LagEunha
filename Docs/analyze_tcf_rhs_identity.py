"""Check frozen-snapshot TCF RHS thermodynamic identity, not time stability."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('root')
ap.add_argument('--output', required=True)
ap.add_argument('--cases', nargs='+', default=['t4', 't6'])
args = ap.parse_args()
root = Path(args.root)
results = {}
for case in args.cases:
    d = root/case
    assert (d/'exit.code').read_text().strip() == '0'
    assert '[TCF] dump written' in (d/'log').read_text()
    p = d/'rhs.dat'
    a = np.loadtxt(p)
    assert a.shape == (36864, 25) and np.isfinite(a).all()
    assert len(np.unique(a[:,0])) == len(a)
    r = np.hypot(a[:,1]-3, a[:,2]-3)
    du = a[:,10]-a[:,3]*a[:,8]-a[:,4]*a[:,9]
    ad = -a[:,5]*a[:,12]
    heat = a[:,11]
    residual = du-ad-heat
    report = {}
    for name, mask in [('inner',r<.5), ('disk',(r>=.5)&(r<=2)),
                       ('outer_active',(r>2)&(r<2.5)), ('all',np.ones(len(a),dtype=bool))]:
        report[name] = dict(n=int(mask.sum()), dU_rate=float(du[mask].sum()),
                            compression_rate=float(ad[mask].sum()),
                            TCF_heat_rate=float(heat[mask].sum()),
                            heat_min=float(heat[mask].min()),
                            residual_abs_sum=float(np.abs(residual[mask]).sum()),
                            residual_max=float(np.abs(residual[mask]).max()))
    results[case] = {'rhs_sha256':hashlib.sha256(p.read_bytes()).hexdigest(), 'regions':report}
result = {'cases':results,
          'scope':'One pre-HG force evaluation per saved state; instantaneous radial regions. Raw force RHS includes frozen cells; all-region values are not actual applied increments. No time-integrated or long-term stability claim.'}
Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

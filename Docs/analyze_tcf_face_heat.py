"""Attribute existing directed-face TCF heat; no evolution changes."""
import argparse
import json
from pathlib import Path
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('root')
ap.add_argument('--output', required=True)
args = ap.parse_args()
root = Path(args.root)
out = {}
for case in ('t4','t6'):
    f = np.loadtxt(root/case/'rhs.dat.faces')
    c = np.loadtxt(root/case/'rhs.dat')
    assert f.shape[1] == 21 and np.isfinite(f).all()
    A, du, Q, Z, qin, delta, chi = (f[:,i] for i in (4,9,14,15,18,19,20))
    mm = np.where(qin*du>0, np.sign(qin)*np.minimum(abs(qin),abs(du)), 0)
    # Current canonical Q law is Z*minmod. Decimal dumps have 7 sig digits.
    assert np.allclose(Q,Z*mm,rtol=3e-5,atol=1e-16)
    assert np.allclose(qin,(1-chi)*delta+chi*du,rtol=3e-5,atol=2e-7)
    heat = .5*A*Q*du
    assert heat.min() >= -1e-20
    r = np.hypot(f[:,2]-3,f[:,3]-3)
    cr = np.hypot(c[:,1]-3,c[:,2]-3)
    out[case] = {}
    for name, m, cm in [('inner',r<.5,cr<.5), ('disk',(r>=.5)&(r<=2),(cr>=.5)&(cr<=2))]:
        total = float(heat[m].sum())
        cell_total = float(c[cm,11].sum())
        assert abs(total-cell_total) < 2e-5*abs(cell_total)
        # Hold impedance and geometry fixed: attribution probe only, not
        # a prediction of a chi-disabled simulation (Z also depends on mm).
        mm0 = np.where(delta*du>0,np.sign(delta)*np.minimum(abs(delta),abs(du)),0)
        out[case][name] = {
            'directed_faces':int(m.sum()), 'heat':total,
            'cell_heat':cell_total,
            'approaching_face_heat_fraction':float(heat[m&(du>0)].sum()/total),
            'receding_face_heat_fraction':float(heat[m&(du<0)].sum()/total),
            'heat_weighted_chi':float((heat[m]*chi[m]).sum()/total),
            'chi_ge_half_heat_fraction':float(heat[m&(chi>=.5)].sum()/total),
            'delta_only_heat_fixed_Z':float((.5*A*Z*mm0*du)[m].sum()),
            'chi_bins':[{'range':[lo,hi], 'heat_fraction':float(heat[m&(chi>=lo)&(chi<hi)].sum()/total)}
                        for lo,hi in [(0,.01),(.01,.1),(.1,.5),(.5,1.000001)]]}
result = {'cases':out, 'scope':'Instantaneous directed-face contributions to owner-cell heat. Approaching faces in shear do not establish physical compression; chi is reconstructed compression share, not an independent truth. Delta-only fixed-Z probe is not an evolved alternative.'}
Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

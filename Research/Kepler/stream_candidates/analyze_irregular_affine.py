from pathlib import Path
import json, sys
import numpy as np
from audit_snapshot import read_snapshot

root = Path(__file__).resolve().parent
base = root / (sys.argv[1] if len(sys.argv)>1 else 'irregular_affine_probe')
manifest = json.loads((base / 'manifest.json').read_text())
results = {}
for case in manifest['cases']:
    d = base / case['name']
    p, meta = read_snapshot(d / 'khout.000002.dat')
    initial, _ = read_snapshot(d / 'khout.000001.dat')
    ids = (p['u4if'] & np.uint64((1 << 56) - 1)).astype(int)
    initial_ids = (initial['u4if'] & np.uint64((1 << 56) - 1)).astype(int)
    lookup = np.empty(len(ids), dtype=int)
    lookup[initial_ids] = np.arange(len(ids))
    f = np.loadtxt(d / 'face_dump_r0_call1.tsv', skiprows=1)
    i = lookup[f[:, 4].astype(int)]
    j = lookup[f[:, 5].astype(int)]
    r = np.hypot(initial['x'] - 12, initial['y'] - 12)
    s, t = manifest['rate'], meta['t']
    if case['flow'] == 'rigid':
        exact = np.array([[s*s*t, -s], [s, s*s*t]]) / (1+s*s*t*t)
    else:
        exact = np.diag([s/(1+s*t), -s/(1-s*t)])
    grad = np.array([[p['stress.gUxx'], p['stress.gUxy']],
                     [p['stress.gUyx'], p['stress.gUyy']]])
    r_out = np.hypot(p['x'] - 12, p['y'] - 12)
    V = initial['volume'].astype(float)
    P = initial['pressure'].astype(float)
    heat_face = f[:,9] + f[:,7] * (initial['vx'][i]*f[:,10] + initial['vy'][i]*f[:,11])
    Ud = np.bincount(i, weights=heat_face, minlength=len(ids))
    Vd = np.bincount(i, weights=f[:,8], minlength=len(ids))
    extra = Ud + P*Vd
    radial = []
    for lo, hi in [(0,.5),(.5,1),(1,2),(2,4),(4,8)]:
        sel = (r >= lo)&(r < hi)
        fs = sel[i]
        gout = (r_out >= lo)&(r_out < hi)
        err = np.sqrt(np.sum((grad[:,:,gout]-exact[:,:,None])**2, axis=(0,1))) / np.linalg.norm(exact)
        radial.append({'Rlo':lo, 'Rhi':hi, 'n':int(sel.sum()),
                       'first_pstar_max_over_P':float(np.max(f[fs,7])/manifest['P']),
                       'first_pstar_rms_rel_error':float(np.sqrt(np.mean((f[fs,7]/manifest['P']-1)**2))),
                       'first_sum_positive_extra_heat_rate':float(np.maximum(extra[sel],0).sum()),
                       'first_volume_rate_over_volume_rms':float(np.sqrt(np.mean((Vd[sel]/V[sel])**2))),
                       'final_grad_relative_error_max':float(err.max()),
                       'final_grad_relative_error_p99':float(np.quantile(err,.99)),
                       'final_grad_relative_error_median':float(np.median(err))})
    results[case['name']] = {'t':t,'flow':case['flow'], 'exact_free_affine_gradient':exact.tolist(),
                             'caution':'first-force metrics diagnose initial operator; final gradient includes one tiny RK step, no centripetal force in this patch test',
                             'radial':radial}
    print(case['name'], json.dumps(radial[0]))
(base / 'results.json').write_text(json.dumps(results,indent=2)+'\n')

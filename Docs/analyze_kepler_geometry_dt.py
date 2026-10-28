"""Reconstruct canonical face timestep limits from existing CPU RHS dumps.

Printed face factors are clamped at 1e-12; omit such faces rather than
pretending their unclamped timestep factors are recoverable.
"""
import argparse
import json
from pathlib import Path
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('run', type=Path)
p.add_argument('--output', required=True, type=Path)
a = p.parse_args()
c = np.loadtxt(a.run / 'rhs.dat')
c = c[np.argsort(c[:, 0])]
f = np.loadtxt(a.run / 'rhs.dat.faces')
ids = c[:, 0].astype(np.int64)
i = np.searchsorted(ids, f[:, 0].astype(np.int64))
j = np.searchsorted(ids, f[:, 1].astype(np.int64))
assert np.all(i < len(c)) and np.all(j < len(c))
assert np.array_equal(ids[i], f[:, 0]) and np.array_equal(ids[j], f[:, 1])
r = np.hypot(c[:, 1]-3, c[:, 2]-3)
active = (r[i] < 2.5) | (r[j] < 2.5)
recoverable = (f[:, 16] > 1e-12) & (f[:, 17] > 1e-12)
gamma = 1.666667
cs = np.sqrt(gamma*c[:, 5]*c[:, 6]/c[:, 7])
sound_error = np.max(np.abs(cs[i]-f[:, 11])/cs[i])
assert sound_error < 2e-6, sound_error
cf = np.maximum(cs[i]*f[:, 16]**((gamma-1)/2),
                cs[j]*f[:, 17]**((gamma-1)/2))
dv = c[j, 3:5]-c[i, 3:5]
normal = np.sum(dv*f[:, 6:8], axis=1)
d = f[:, 5]
dtg = .125*d/(cf+np.abs(normal))
cfl = .8*d/(cs[i]+cs[j]-np.minimum(0, normal))
valid = np.flatnonzero(active & recoverable)
order = valid[np.argsort(dtg[valid])[:10]]
rows = []
for k in order:
    rows.append(dict(i=int(f[k, 0]), j=int(f[k, 1]), ri=float(r[i[k]]),
        rj=float(r[j[k]]), separation=float(d[k]), dt_geom=float(dtg[k]),
        dt_cfl=float(cfl[k]), cs_i=float(cs[i[k]]), cs_j=float(cs[j[k]]),
        face_sound_max=float(cf[k]), normal_relative_speed=float(abs(normal[k])),
        sound_denominator_fraction=float(cf[k]/(cf[k]+abs(normal[k]))),
        fi=float(f[k, 16]), fj=float(f[k, 17]),
        separation_over_sqrt_volume_i=float(d[k]/np.sqrt(c[i[k], 6]))))
out = dict(source=str(a.run), gamma=gamma,
    scope='instantaneous pre-HG RHS; approximate reconstruction from printed precision',
    limitation='Only emitted TCF faces are represented; this is not a reconstruction of every limiter evaluation or subsequent RK stage.',
    active_directed_faces=int(active.sum()),
    excluded_clamped_active_faces=int((active & ~recoverable).sum()),
    reconstructed_cell_sound_max_relative_error=float(sound_error),
    shortest_active_separation=float(d[active].min()), limiting_faces=rows)
with a.output.open('x') as stream:
    json.dump(out, stream, indent=2)
    stream.write('\n')
print(json.dumps(out, indent=2))

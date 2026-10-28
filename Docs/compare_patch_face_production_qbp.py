"""Read-only saved t2 cache comparison; not a production HLLC trace."""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read

base = Path('/gpfs/kjhan/LagForce')
root = base / 'Codex_pair_hg_modes_t2_20261018'
manifest = json.loads((root / 'manifest.json').read_text())
source = Path(manifest['source'])
assert hashlib.sha256(source.read_bytes()).hexdigest() == manifest['source_sha256']
ic, time = read(source)
failure = base / 'Codex_checkpoint_patch_rhs_exact_failure_r0.5_20261018'
frames = json.loads((failure / 'report.json').read_text())['results'][0]['failure_frames']
a, b = frames[0]['values'], frames[1]['values']
saved = np.load(failure / 'source.npz')
ids = [int(saved['ids'][a[k]]) for k in ('i', 'j')]
cache = np.loadtxt(root / 'raw_faces.dat.cache.rank0000')
rows = np.array([cache[cache[:, 0] == k][0] for k in ids])
assert np.all(rows[:, 2] == 1) and np.all(rows[:, 5] != 0)
x = np.array([[ic['x'][k]-3, ic['y'][k]-3] for k in ids])
patch_positions = saved['state'][:-4].reshape(-1, 20)[:, :2]
assert np.array_equal(x, patch_positions[[a['i'], a['j']]])
v = np.array([[ic['vx'][k], ic['vy'][k]] for k in ids])
G = rows[:, 7:11].reshape(2, 2, 2)
H = rows[:, 11:17].reshape(2, 2, 3)
T = rows[:, 17:25].reshape(2, 2, 4)
n = np.array(a['n']); point = np.array(a['point'])
traces = {}
for order in (1, 2, 3):
    out = v.copy()
    for k in range(2):
        rx, ry = point-x[k]
        out[k] += G[k] @ np.array([rx, ry])
        if order >= 2:
            out[k] += .5*(H[k] @ np.array([rx*rx, 2*rx*ry, ry*ry]))
        if order >= 3:
            out[k] += T[k] @ np.array([rx**3, 3*rx*rx*ry, 3*rx*ry*ry, ry**3])/6
    traces[str(order)] = float((out[1]-out[0]) @ n)
d = np.linalg.norm(x[1]-x[0]); e = (x[1]-x[0])/d
arc = rows[:, 6].astype(bool)
effective = G.copy()
if arc[0] != arc[1]:
    inner = 1 if arc[0] else 0
    rx, ry = (-d if arc[0] else d)*e
    ext = G[inner].copy()
    for c in range(2):
        h, t = H[inner, c], T[inner, c]
        ext[c, 0] += h[0]*rx+h[1]*ry+.5*(t[0]*rx*rx+2*t[1]*rx*ry+t[2]*ry*ry)
        ext[c, 1] += h[1]*rx+h[2]*ry+.5*(t[1]*rx*rx+2*t[2]*rx*ry+t[3]*ry*ry)
    effective[1-inner] = ext
ducell = float((v[0]-v[1]) @ e)
qdel = float(ducell+.5*d*(e @ effective.sum(axis=0) @ e))
chis = []
for g in effective:
    norm = np.linalg.norm(.5*(g+g.T))
    chis.append(min(1., max(-np.trace(g), 0.)/norm) if norm else 0.)
chi = max(chis); qin = (1-chi)*qdel+chi*ducell
mm = 0. if qin*ducell <= 0 else min((qin, ducell), key=abs)
faces = np.loadtxt(root / 'raw_faces.dat.rank0000')
match = ((faces[:, 0] == ids[0]) & (faces[:, 1] == ids[1])) | ((faces[:, 0] == ids[1]) & (faces[:, 1] == ids[0]))
assert np.any(match)
print(json.dumps(dict(source=manifest, actual_time=float(time), ids=ids,
    original_face_rows=faces[match].tolist(), arc=arc.tolist(),
    normal_alignment=float(e@n), point_normal_jump=float((v[1]-v[0])@n),
    conditional_affine_jump=b['ur']-b['ul'],
    conditional_vacuum_threshold=2*(b['cl']+b['cr'])/(b['gamma']-1),
    unlimited_polynomial_normal_jumps=traces,
    production_qbp_formula=dict(ducell=ducell, qdel=qdel, chi=chi, qin=qin, minmod=mm),
    scope='Saved production cache; QBP EDGE_INNER=3/CHI=1 formula replay. Polynomial evaluations are unlimited diagnostics at remapped face, not production Riemann traces. QBP minmod is an impedance input, not a face velocity jump. No evolution or physical vacuum claim.'), indent=2))

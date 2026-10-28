"""Diagnostic uniform-density affine lift; no state modification or evolution.
Hold original cell mass, mean velocity and hydro energy; charge gradient kinetic
energy against stored ie. This is a feasibility check, not an approved restart.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

base = Path('/gpfs/kjhan/LagForce')
sys.path.insert(0, str(base))
sys.path.insert(0, str(base/'Codex_coupled_relaxation_evolution_20261018'))
from run_canonical_av_branches_20261018 import read
from periodic2d_hllc import quadrature

root = base/'Codex_pair_hg_modes_t2_20261018'
manifest = json.loads((root/'manifest.json').read_text())
source = Path(manifest['source'])
assert hashlib.sha256(source.read_bytes()).hexdigest() == manifest['source_sha256']
ic, time = read(source)
x = np.column_stack([ic['x']-3, ic['y']-3])
faces = np.loadtxt(root/'raw_faces.dat.rank0000')
cache = np.loadtxt(root/'raw_faces.dat.cache.rank0000')
cache = cache[np.argsort(cache[:, 0])]
assert np.array_equal(cache[:, 0], np.arange(len(ic)))
order = np.argsort(faces[:, 0], kind='stable'); owners = faces[order, 0]
ids = np.flatnonzero(np.linalg.norm(x, axis=1)<.7)
records = []
for i in ids:
    rows = faces[order[np.searchsorted(owners, i):np.searchsorted(owners, i, side='right')]]
    assert len(rows) >= 3 and cache[i, 5] != 0
    vertices = rows[:, 6:8]  # generator-relative, avoids subtracting large origin
    angles = np.arctan2(vertices[:, 1], vertices[:, 0])
    pts, dv = quadrature(vertices[np.argsort(angles)])
    volume = dv.sum(); center = dv@pts/volume
    dm = dv*(ic['mass'][i]/volume)
    z = pts-center
    S = z.T@(dm[:, None]*z)
    G = cache[i, 7:11].reshape(2, 2)
    K = .5*np.trace(G@S@G.T)
    C = S@G.T; ell = C[0, 1]-C[1, 0]
    Krot = ell*ell/(2*np.trace(S))
    # Moving the mass centroid from generator to polygon centroid changes
    # orbital L even with the same mean momentum; report this separately.
    p = ic['mass'][i]*np.array([ic['vx'][i], ic['vy'][i]])
    orbital_change = center[0]*p[1]-center[1]*p[0]
    U = float(ic['ie'][i])
    assert U > 0 and K >= 0 and Krot <= K+1e-18
    records.append(dict(id=int(i), radius=float(np.linalg.norm(x[i])),
        volume_relative_error=float(volume/ic['volume'][i]-1),
        kinetic=float(K), rotational=float(Krot), thermal=U,
        ratio=float(K/U), remaining_thermal=float(U-K),
        angular_change=float(orbital_change+ell)))
out = {}
for name, select in [('all_r_lt_07', lambda r: True),
                     ('inner_r_lt_05', lambda r: r['radius']<.5),
                     ('interface_045_055', lambda r: .45<r['radius']<.55)]:
    rr = [r for r in records if select(r)]
    ratios = np.array([r['ratio'] for r in rr])
    out[name] = dict(cells=len(rr), nonpositive_thermal=sum(r['remaining_thermal']<=0 for r in rr),
        ratio_quantiles=np.quantile(ratios, [0, .5, .9, .99, 1]).tolist(),
        sum_kinetic=sum(r['kinetic'] for r in rr), sum_thermal=sum(r['thermal'] for r in rr),
        sum_angular_change=sum(r['angular_change'] for r in rr),
        maximum_volume_relative_error=max(abs(r['volume_relative_error']) for r in rr))
print(json.dumps(dict(scope=__doc__, manifest=manifest, time=float(time), regions=out,
    failed_face=[r for r in records if r['id'] in (755,864)],
    worst_cells=sorted(records, key=lambda r:r['ratio'], reverse=True)[:10],
    limitations='Uniform density on original unweighted polygons; gradient centered to preserve mean p. Does not preserve original angular momentum automatically, pressure, or generator velocity. Stored ie used directly including any original subcell convention. No general impossibility claim for alternative coupled representations.'), indent=2))

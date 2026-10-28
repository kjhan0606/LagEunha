"""Kepler density fields from stored cell means, without particle markers.

Unweighted Voronoi geometry is valid only when every stored weight is zero.
No spatial smoothing is applied. Run on a compute host.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
import numpy as np
from scipy.spatial import Voronoi

base = Path('/gpfs/kjhan/LagForce')
sys.path.insert(0, str(base))
from run_canonical_av_branches_20261018 import read

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--snapshot', type=Path, required=True)
parser.add_argument('--output-prefix', type=Path, required=True)
args = parser.parse_args()
png = Path(str(args.output_prefix)+'.png')
meta_path = Path(str(args.output_prefix)+'.json')
assert not png.exists() and not meta_path.exists()
icpath = base/'Codex_canonical_initial_20261017/run/khout.000000.dat'
ic, t0 = read(icpath)
data, t = read(args.snapshot)
assert t0 == 0 and np.all(ic['w2'] == 0) and np.all(data['w2'] == 0)
assert np.array_equal(ic['mass'], data['mass'])
assert np.array_equal(ic['id'] & np.uint64(2**56-1), data['id'] & np.uint64(2**56-1))
fig, axs = plt.subplots(2, 2, figsize=(12, 11), layout='constrained')
norm = LogNorm(.005, 2.)
counts = []
for col, (frame, time) in enumerate(((ic, t0), (data, t))):
    x = np.column_stack([frame['x']-3, frame['y']-3])
    assert np.isfinite(x).all() and np.all(np.isfinite(frame['den'])) and np.all(frame['den'] > 0)
    vor = Voronoi(x)
    polys, rho = [], []
    for i, index in enumerate(vor.point_region):
        region = vor.regions[index]
        if not region or -1 in region:
            assert not np.all(abs(x[i]) < 2.4), 'Unbounded cell near visible domain'
            continue
        polys.append(vor.vertices[region])
        rho.append(frame['den'][i])
    counts.append(len(polys))
    for row in range(2):
        ax = axs[row, col]
        pc = PolyCollection(polys, array=np.asarray(rho), cmap='viridis', norm=norm,
                            edgecolors='none', linewidth=0, antialiaseds=False, rasterized=True)
        ax.add_collection(pc)
        limit = 2.2 if row == 0 else .65
        ax.set(xlim=(-limit, limit), ylim=(-limit, limit), xlabel='x - center', ylabel='y - center')
        ax.set_aspect('equal')
        ax.set_title(('Initial' if col == 0 else f't = {time:.5f}')+
                     (' density field' if row == 0 else ' density field: inner zoom'))
fig.colorbar(pc, ax=axs, label='Cell-averaged density', shrink=.86)
fig.suptitle('Kepler density field: existing corotating auxiliary-memory run B\n'
             'Stored cell-average density on current unweighted Voronoi cells', fontsize=13)
png.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(png, dpi=180)
plt.close(fig)
meta = dict(scope=__doc__, snapshot=str(args.snapshot), time=t, image=str(png),
    source_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__), args.snapshot, icpath)},
    finite_polygon_counts=counts, density_color_range=[.005, 2.], weights_zero=True,
    rendering='Piecewise constant stored cell density, no point markers, cell outlines or interpolation. Unbounded outer cells outside displayed domain omitted. Not new moment model.')
with meta_path.open('x') as stream:
    json.dump(meta, stream, indent=2)
print(json.dumps(meta, indent=2))

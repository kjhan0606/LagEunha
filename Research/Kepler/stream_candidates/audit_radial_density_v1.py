"""Actual W0 density fields and independent cell-average reference at t=.05."""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

from audit_snapshot import read_snapshot
from prepare_radial_resolution_v1 import geometry, resolve_near_origin_edges, integrate_density

parser = argparse.ArgumentParser()
parser.add_argument('--root', type=Path, required=True)
args = parser.parse_args()
base = args.root/'radial_resolution_20261021_v2'
result = dict(scope='matched t=.05 W0 native checkpoints; independent stationary density integrated over CURRENT true cells',
              long_goal_achieved=False, cases={})
fig, axes = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)
for column, label in enumerate(('uniform', 'graded')):
    native = json.loads((base/'native_v2'/label/'manifest.json').read_text())
    stage = next(x for x in native['stages'] if x['label'] == 't0.05' and x['exit_code'] == 0)
    path = Path(stage['checkpoint'])
    p, meta = read_snapshot(path)
    if meta['t'] != .05 or np.any(p['w2'] != 0):
        raise ValueError('requires matched t=.05 W0 snapshots')
    xy = np.column_stack((p['x'], p['y']))
    fan, _ = geometry(xy)
    true_volume = fan[0]
    fan = resolve_near_origin_edges(fan)
    coarse = integrate_density(4, fan)
    reference = integrate_density(8, fan)
    error = abs(reference-coarse)/reference
    for order in (16, 32, 64):
        bad = error > 1e-8
        if not bad.any():
            break
        revised = integrate_density(order, fan, bad)
        error[bad] = abs(revised[bad]-reference[bad])/revised[bad]
        reference[bad] = revised[bad]
    if error.max() > 1e-8:
        raise ValueError('current-cell density reference quadrature incomplete')
    density = p['mass']/p['volume']
    relative = abs(p['mass']-reference)/reference
    radius = np.linalg.norm(xy-12, axis=1)
    row = dict(checkpoint=str(path), time=meta['t'], n=len(p),
               native_volume_relative_max=float(np.max(abs(p['volume']-true_volume)/true_volume)),
               reference_quadrature_relative_change_max=float(error.max()), radial=[])
    for lo, hi in ((0, .047), (.047, .1), (.1, .2), (.2, .5), (.5, 1), (1, 2), (2, 4), (4, 8)):
        sel = (radius >= lo) & (radius < hi)
        if not sel.any():
            continue
        row['radial'].append(dict(current_Rlo=lo, current_Rhi=hi, n=int(sel.sum()),
            density_relative_L1_mass_weighted=float(np.sum(p['mass'][sel]*relative[sel])/p['mass'][sel].sum()),
            density_relative_max=float(relative[sel].max()),
            density_min=float(density[sel].min()), density_max=float(density[sel].max())))
    result['cases'][label] = row
    tree = cKDTree(xy % 24, boxsize=24)
    for panel, extent in enumerate((9., .5)):
        axis = axes[panel, column]
        coordinate = (np.arange(768)+.5)/768*2*extent-extent
        xx, yy = np.meshgrid(coordinate, coordinate)
        query = np.column_stack((xx.ravel()+12, yy.ravel()+12))
        indices = tree.query(query, workers=1)[1]
        image = density[indices].reshape(xx.shape)
        norm = LogNorm(1e-3, 1) if panel == 0 else Normalize(.0008, .0012)
        artist = axis.imshow(image, origin='lower', extent=(-extent, extent, -extent, extent),
                             cmap='viridis', norm=norm, interpolation='nearest')
        axis.set_title(f'{label}, N={len(p):,}, t=0.05'+(' (core)' if panel else ''))
        axis.set_xlabel('x - center')
        axis.set_ylabel('y - center')
        fig.colorbar(artist, ax=axis, label='surface density', extend='both')
result['visualization'] = dict(field='piecewise constant mass/stored_volume in actual unweighted Voronoi cells',
    sampling='periodic nearest-generator queries; no scatter interpolation',
    core_color_limits=[.0008, .0012], disk_color_limits=[.001, 1],
    clipped_color_values='colorbars extend at both ends; quantitative extrema are in the JSON')
(base/'density_t005.json').write_text(json.dumps(result, indent=2)+'\n')
fig.savefig(base/'density_t005.png', dpi=180)
plt.close(fig)
print(json.dumps(result, indent=2))

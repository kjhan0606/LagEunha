"""Compare bounded Q1 saved density fields at shared fixed color limits."""
import argparse
import json
from pathlib import Path
import socket
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from plot_kepler_q1_boundary_comparison_20261029 import sample, radial, digest


def plot(initial, endpoints, output, summary):
    if socket.gethostname().split('.')[0].lower() in ('grammar', 'syntax'):
        raise RuntimeError('Render on a compute host')
    if output.exists() or summary.exists():
        raise FileExistsError('Refusing to overwrite a density history')
    a, m0 = ordered_global(initial); ids = a['u4if'] & ID_MASK
    if m0['t'] != 0: raise ValueError('Requires exact original t0 reference')
    r0 = np.hypot(a['x']-3., a['y']-3.); disk = (r0 >= .5) & (r0 < 2.)
    states = [(initial, a, m0)]
    for endpoint in endpoints:
        b, mt = ordered_global(endpoint, ids)
        if mt['t'] <= states[-1][2]['t'] or not np.array_equal(a['mass'], b['mass']):
            raise ValueError('Require increasing times and fixed mass/IDs')
        states.append((endpoint, b, mt))
    fig, axes = plt.subplots(2, len(states), figsize=(5*len(states), 9),
                             squeeze=False, constrained_layout=True)
    rows = []
    for column, (path, p, meta) in enumerate(states):
        if np.any(p['w2'] != 0) or not np.isfinite(p['den']).all() or np.any(p['den'] <= 0):
            raise ValueError('Requires finite positive W0 density')
        delta, profile = radial(a, p)
        images = [sample(p, p['den'], bounded=True),
                  sample(p, 100*delta, bounded=True, disk_mask=disk)]
        for row, values in enumerate(images):
            ax = axes[row, column]
            im = ax.imshow(values, origin='lower', extent=(-2.25,2.25,-2.25,2.25),
                           norm=Normalize(0,1.05) if row == 0 else Normalize(-1,1),
                           cmap='viridis' if row == 0 else 'coolwarm')
            ax.set(title=f't = {meta["t"]:g}', xlabel='x - center', ylabel='y - center', aspect='equal')
        rows.append({'time': meta['t'], 'source': str(path), 'sha256': digest(path),
                     'radial_disk': profile})
    for row, label in enumerate(['Surface density', 'Initial-ID disk density change (%)']):
        fig.colorbar(axes[row, -1].images[0], ax=axes[row, :].tolist(), label=label, extend='both', shrink=.9)
    fig.suptitle('Native bounded Q1, W=0: common scales, piecewise-constant cell density')
    output.parent.mkdir(parents=True, exist_ok=True); fig.savefig(output, dpi=170); plt.close(fig)
    result = {'initial_sha256': digest(initial), 'density_limits': [0,1.05],
              'density_change_percent_limits': [-1,1], 'disk_cells': int(disk.sum()), 'states': rows,
              'scope': 'W0 nearest-site cell-average density without smoothing, clipped to the polygon. Lower panels display only initial disk IDs; white regions exclude other IDs and do not measure vacuum. Colors saturate outside fixed limits. Images are morphology evidence, not a substitute for physical gates or GIZMO long-time validation.'}
    summary.write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial', type=Path, required=True)
    parser.add_argument('--endpoints', type=Path, nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--summary', type=Path, required=True)
    print(json.dumps(plot(**vars(parser.parse_args())), indent=2))

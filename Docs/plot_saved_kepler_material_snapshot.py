"""Render real saved material-point Kepler data on a compute host.

Stored cell values color generator points; this is not polygon rendering.
Material phase diagnostics apply to the old auxiliary-memory model only.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

base = Path('/gpfs/kjhan/LagForce')
sys.path.insert(0, str(base))
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--snapshot', type=Path, required=True)
parser.add_argument('--output-prefix', type=Path, required=True)
parser.add_argument('--label', default='Existing corotating auxiliary-memory run B')
args = parser.parse_args()
png = Path(str(args.output_prefix)+'.png')
meta_path = Path(str(args.output_prefix)+'.json')
assert not png.exists() and not meta_path.exists(), 'Preserve existing figures'
icpath = base/'Codex_canonical_initial_20261017/run/khout.000000.dat'
ic, t0 = read(icpath)
a, t = read(args.snapshot)
assert t0 == 0. and np.isfinite(t)
assert np.array_equal(a['mass'], ic['mass'])
assert np.array_equal(a['id'] & np.uint64(2**56-1), ic['id'] & np.uint64(2**56-1))
x0 = np.column_stack([ic['x']-3, ic['y']-3])
x = np.column_stack([a['x']-3, a['y']-3])
r0 = np.linalg.norm(x0, axis=1)
phase = np.angle(np.exp(1j*(np.arctan2(x[:, 1], x[:, 0])-np.arctan2(x0[:, 1], x0[:, 0])
                          -np.sqrt(gravity_coefficient(r0))*t)))*180/np.pi
perr = np.abs(a['pressure']/ic['pressure']-1)
assert all(np.isfinite(v).all() for v in (x0, x, phase, perr, a['den']))
assert np.all(a['den'] > 0)
fig, axs = plt.subplots(2, 2, figsize=(12, 11), layout='constrained')
for ax, pos, rho, title in ((axs[0, 0], x0, ic['den'], 'Initial density, t = 0'),
                            (axs[0, 1], x, a['den'], f'Density, t = {t:.5f}')):
    scatter = ax.scatter(pos[:, 0], pos[:, 1], c=rho, s=2.3, cmap='viridis',
                         norm=LogNorm(.01, 2), rasterized=True, linewidths=0)
    ax.set_title(title)
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-2.2, 2.2)
    fig.colorbar(scatter, ax=ax, label='Stored cell density')
inner = np.linalg.norm(x, axis=1) < .75
ax = axs[1, 0]
scatter = ax.scatter(x[inner, 0], x[inner, 1], c=np.maximum(perr[inner], 1e-8), s=9,
    cmap='magma', norm=LogNorm(1e-5, 30), rasterized=True, linewidths=0)
ax.set_title('Inner region: absolute relative pressure error')
fig.colorbar(scatter, ax=ax, label='|P / P(initial) - 1|')
ax = axs[1, 1]
scatter = ax.scatter(x[inner, 0], x[inner, 1], c=phase[inner], s=9, cmap='coolwarm',
                     vmin=-180, vmax=180, rasterized=True, linewidths=0)
ax.set_title('Inner region: wrapped material-orbit phase error')
fig.colorbar(scatter, ax=ax, label='Phase error [degrees]')
for ax in axs.flat:
    ax.set_aspect('equal')
    ax.set_xlabel('x - center')
    ax.set_ylabel('y - center')
    for radius in (.5, 2.):
        ax.add_patch(plt.Circle((0, 0), radius, fill=False, edgecolor='gray', ls='--', lw=.7, alpha=.8))
for ax in axs[1]:
    ax.set_xlim(-.7, .7)
    ax.set_ylim(-.7, .7)
fig.suptitle(f'Kepler disk: {args.label}\nActual saved data at t = {t:.8f}; not the new moment model', fontsize=13)
png.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(png, dpi=170)
plt.close(fig)
meta = dict(scope=__doc__, snapshot=str(args.snapshot), initial=str(icpath), time=t, image=str(png),
    source_sha256={str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in (Path(__file__), args.snapshot, icpath)},
    color_ranges=dict(density=[.01, 2.], absolute_relative_pressure_error=[1e-5, 30], phase_degrees=[-180, 180]),
    rendering='Generator point scatter; color scales clip outside displayed ranges. No synthesized data or new evolution.',
    numpy=np.__version__, matplotlib=matplotlib.__version__)
with meta_path.open('x') as stream:
    json.dump(meta, stream, indent=2)
print(json.dumps(meta, indent=2))

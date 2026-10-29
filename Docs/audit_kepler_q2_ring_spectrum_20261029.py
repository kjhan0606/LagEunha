"""Angular spectrum of specific angular-momentum residual on original Q2 rings.

Uses existing bounded-domain snapshots; Fourier periodicity is solely the
cyclic ordering of material IDs around a ring, not a spatial boundary mode.
"""
import hashlib
import json
import socket
from pathlib import Path

import numpy as np
from kepler_native_global_ids_20261029 import ID_MASK, ordered_global

if socket.gethostname().split('.')[0].lower() in ('grammar', 'syntax'):
    raise RuntimeError('Run analysis on a compute host')

docs = Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
late = json.loads((docs / 'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
feedback = json.loads((docs / 'KEPLER_Q2_TORQUE_FEEDBACK_2026-10-29.json').read_text())
root = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()

initial = Path(late['initial'])
a, meta = ordered_global(initial)
if meta['t'] != 0 or digest(initial) != late['initial_sha256']:
    raise ValueError('Initial snapshot differs')
ids = a['u4if'] & ID_MASK
x = a['x'].astype(float) - 3
y = a['y'].astype(float) - 3
r = np.hypot(x, y)
angle = np.arctan2(y, x)
L0 = x * a['vy'].astype(float) - y * a['vx'].astype(float)
selected = (r >= 0.5) & (r < 0.75)
rings = []
for radius in np.unique(np.round(r[selected], 10)):
    ix = np.flatnonzero(selected & (np.round(r, 10) == radius))
    if len(ix) < 24:
        raise ValueError('Unexpected incomplete ring')
    ix = ix[np.argsort(angle[ix])]
    gaps = np.diff(np.r_[angle[ix], angle[ix][0] + 2 * np.pi])
    if np.max(np.abs(gaps / (2 * np.pi / len(ix)) - 1)) > 1e-4:
        raise ValueError('Ring not uniformly spaced')
    rings.append(ix)

paths = [Path(late['cases'][0]['endpoint'])] + [root / ('khout.' + n + '.dat') for n in
         ('001329', '002762', '004411', '005767', '007138', '008700', '011058')]
cases = []
for path in paths:
    p, meta = ordered_global(path, ids)
    sha = digest(path)
    reference = next((item for item in feedback['cases'] if item['time'] == meta['t']), None)
    if reference is None or sha != reference['source_sha256'] or not np.array_equal(a['mass'], p['mass']) or np.any(p['w2'] != 0):
        raise ValueError('Unexpected snapshot')
    xp = p['x'].astype(float) - 3
    yp = p['y'].astype(float) - 3
    residual = xp * p['vy'].astype(float) - yp * p['vx'].astype(float) - L0
    power = np.zeros(5)
    peak_modes = []
    for ix in rings:
        values = residual[ix] - np.mean(residual[ix])
        ft = np.fft.rfft(values) / len(ix)
        spectrum = np.abs(ft)**2
        spectrum[1:] *= 2
        if len(ix) % 2 == 0:
            spectrum[-1] *= 0.5
        freq = np.fft.rfftfreq(len(ix))
        # Sum ring power weighted by its site count (all Q2 ring sites have
        # equal mass only within each ring); avoids weighting each ring equally.
        mass = float(np.sum(a['mass'][ix].astype(float)))
        band = np.minimum((freq * 10).astype(int), 4)
        for k in range(5):
            power[k] += mass * np.sum(spectrum[band == k])
        peak_modes.append({'n': len(ix), 'peak_m': int(np.argmax(spectrum[1:]) + 1)})
    if not np.sum(power) > 0:
        raise ValueError('Zero residual power')
    cases.append({'time': meta['t'], 'source_sha256': sha,
                  'band_fractions': (power / np.sum(power)).tolist(),
                  'per_ring_peak_modes': peak_modes})

result = {'initial_sha256': late['initial_sha256'], 'ring_count': len(rings),
          'site_count': int(np.sum(selected)), 'bands_m_over_n': ['[0,.1)', '[.1,.2)', '[.2,.3)', '[.3,.4)', '[.4,.5]'],
          'cases': cases, 'script_sha256': digest(Path(__file__)),
          'scope': 'All original inner-disk rings, mass-weighted power of per-ring mean-subtracted specific angular momentum change from initial state. Cyclic DFT is along ring indices, independent of physical outer boundary. No claim of causal linear growth, temporal mode tracking, or universal wavelength.'}
out = docs / 'KEPLER_Q2_RING_SPECTRUM_2026-10-29.json'
with out.open('x') as stream:
    json.dump(result, stream, indent=2)
    stream.write('\n')
print(json.dumps(result, indent=2))

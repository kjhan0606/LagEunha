"""Prepare native affine patch tests on saved geometry; no fluid evolution here."""
from pathlib import Path
import hashlib, json, struct
import numpy as np
from audit_snapshot import read_snapshot

root = Path(__file__).resolve().parent
out = root / 'irregular_affine_probe'
out.mkdir(exist_ok=False)
source = root / 'final_gpu_late/khout.009571.dat'
template = root / 'affine_probe/khout.000001.dat'
binary = root / 'final_volume_gpu/eunha2_final_volume.exe'
q, meta = read_snapshot(source)
t, _ = read_snapshot(template)
params = (root / 'affine_gpu_lsf/params.dat').read_text()
gamma, pressure, rate = 1.666667, 1.e-6, 0.1
cases = [('rigid_regular_lsf', t, 'rigid', 1, '.333333'),
         ('rigid_irregular_lsf', q, 'rigid', 1, '.333333'),
         ('rigid_irregular_gg0', q, 'rigid', 0, '0.000000'),
         ('rigid_irregular_gg', q, 'rigid', 0, '.333333'),
         ('strain_irregular_lsf', q, 'strain', 1, '.333333')]
manifest = {'purpose': 'native instantaneous affine reconstruction on regular and actual distorted cells; not Kepler stability or intrinsic spin evolution',
            'geometry_snapshot': str(source), 'geometry_time': meta['t'],
            'geometry_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
            'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
            'P': pressure, 'rho': 1, 'rate': rate,
            'gravity': 'disabled by EUNHA_IC=affine_probe',
            'dt': 'small Courant, one native RK4 step; only first-force face dump is initial operator',
            'cases': []}
for label, geometry, flow, method, accuracy in cases:
    dest = out / label
    dest.mkdir()
    raw = np.zeros(len(geometry) * meta['psize'], dtype=np.uint8)
    p = np.ndarray(len(geometry), dtype=q.dtype, buffer=raw)
    for field in ('u4if', 'x', 'y', 'volume'):
        p[field] = geometry[field]
    p['mass'] = geometry['volume']
    p['den'] = 1.
    p['pressure'] = pressure
    p['ie'] = pressure * geometry['volume'] / (gamma - 1)
    dx, dy = p['x'] - 12, p['y'] - 12
    p['vx'] = -rate * dy if flow == 'rigid' else rate * dx
    p['vy'] = rate * dx if flow == 'rigid' else -rate * dy
    fp = dest / 'khout.000001.dat'
    with fp.open('wb') as f:
        f.write(struct.pack('<4i', -1, len(p), meta['psize'], 1))
        f.write(raw.tobytes())
        # KH trailer is t, dt, GAS_dtold, nx, ny, Lx, Ly, alpha, beta.
        # Gamma belongs to params.dat, not to the previous-timestep slot.
        f.write(struct.pack('<ddfiiddff', 0., 1.e-7, 1.e-7, 256, 256, 24., 24., 0., 0.))
    lines = []
    for line in params.splitlines():
        if line.startswith('define GAS gradient_method'):
            line = f'define GAS gradient_method= {method}'
        elif line.startswith('define KH Order of Accuracy'):
            line = f'define KH Order of Accuracy= {accuracy}'
        lines.append(line)
    (dest / 'params.dat').write_text('\n'.join(lines) + '\n')
    (dest / 'eunha2.exe').symlink_to(binary)
    (dest / 'OWNER.txt').write_text('Codex native affine operator diagnostic, 2026-10-20\n')
    manifest['cases'].append({'name': label, 'flow': flow, 'gradient_method': method,
                              'OoA': float(accuracy), 'initial_sha256': hashlib.sha256(fp.read_bytes()).hexdigest()})
(out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps(manifest, indent=2))

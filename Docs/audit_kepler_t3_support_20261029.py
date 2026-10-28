"""Saved-state geometry/support attribution; never evolves or edits a trajectory."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.spatial import Voronoi

from audit_kepler_native_radial_20261026 import ordered, mean, GAMMA_MINUS_ONE


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def geometry(xy, selected, foreign, periodic):
    """For periodic data certify that no omitted image can cut any disk cell.

    All generators must first be inside the box. For a polygon inside a
    radius R ball about its generator, any omitted generator more than 2R
    away cannot cut it. Distance from the generator to the box complement
    bounds the distance to all omitted images. Isolated counterfactuals do
    not use this certificate or claim periodic equivalence.
    """
    if not np.isfinite(xy).all():
        raise ValueError('Nonfinite generator coordinates')
    if periodic and not np.all((xy >= 0.) & (xy <= 6.)):
        raise ValueError('Periodic certificate requires every generator inside the box')
    vor = Voronoi(xy, qhull_options='Qbb Qc Qz')
    areas = np.empty(len(selected))
    radii = np.empty(len(selected))
    polygons = []
    for k, i in enumerate(selected):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError('Unbounded selected cell')
        vertices = vor.vertices[region]
        # Translate before shoelace products to avoid large-coordinate cancellation.
        v = vertices-xy[i]
        areas[k] = .5*abs(np.dot(v[:, 0], np.roll(v[:, 1], -1))-
                         np.dot(v[:, 1], np.roll(v[:, 0], -1)))
        radii[k] = np.max(np.linalg.norm(v, axis=1))
        polygons.append(vertices)
    if not np.all(np.isfinite(areas) & (areas > 0) & (radii > 0)):
        raise ValueError('Invalid independent polygon')
    safety = None
    if periodic:
        edge = np.min(np.column_stack((xy[selected], 6.-xy[selected])), axis=1)
        safety = float(np.min(edge/(2*radii)))
        if safety <= 1.:
            raise ValueError('Omitted periodic images not excluded; require explicit periodic reconstruction')
    local_index = np.full(len(xy), -1, dtype=int)
    local_index[selected] = np.arange(len(selected))
    has_foreign = np.zeros(len(selected), dtype=bool)
    pairs = vor.ridge_points
    support_rows = []
    for a, b in [(pairs[:, 0], pairs[:, 1]), (pairs[:, 1], pairs[:, 0])]:
        keep = (local_index[a] >= 0) & foreign[b]
        owned, other = a[keep], b[keep]
        has_foreign[local_index[owned]] = True
        support_rows.extend(zip(owned.tolist(), other.tolist()))
    return areas, has_foreign, support_rows, safety, polygons


def orbit(xy0, t):
    v = xy0-3.
    r = np.hypot(v[:, 0], v[:, 1])
    phase = t*(r*r+.01175**2)**(-.75)
    return 3.+np.column_stack((v[:, 0]*np.cos(phase)-v[:, 1]*np.sin(phase),
                              v[:, 0]*np.sin(phase)+v[:, 1]*np.cos(phase)))


def row(mask, masses, values):
    if not np.any(mask):
        return {'cells': 0}
    return {'cells': int(mask.sum()), **{key: mean(value[mask], masses[mask])
            for key, value in values.items()}}


def render(path, polygons0, polygons, frozen_polygons, rho0, rho, rho_frozen):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize, LogNorm
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.5), constrained_layout=True)
    for col, (title, vertices, density) in enumerate([
            ('Initial disk cells', polygons0, rho0),
            ('Native disk at t=3', polygons, rho),
            ('Geometry only: exterior restored to t=0', frozen_polygons, rho_frozen)]):
        top = PolyCollection(vertices, array=density, edgecolors='none',
                             cmap='viridis', norm=LogNorm(.01, 1.1))
        axes[0, col].add_collection(top)
        axes[0, col].set_title(title, fontsize=10)
        bottom = PolyCollection(vertices, array=density/rho0, edgecolors='none',
                                cmap='coolwarm', norm=Normalize(.7, 1.3))
        axes[1, col].add_collection(bottom)
        for panel in axes[:, col]:
            panel.set(xlim=(.9, 5.1), ylim=(.9, 5.1), aspect='equal', xlabel='x', ylabel='y')
    fig.colorbar(top, ax=axes[0, :], label='Surface density; shared scale')
    fig.colorbar(bottom, ax=axes[1, :], label='Cell density / its initial density; shared scale')
    fig.suptitle('Saved 384² native CUDA+MPI+OpenMP disk: t=3 radial gate failed\n'
                 'Exterior counterfactual is diagnostic geometry, not evolved fluid')
    fig.savefig(path, dpi=160)
    plt.close(fig)


def analyze(initial, checkpoints, output):
    started = time.monotonic()
    p0, meta0 = ordered(initial)
    if meta0['t'] != 0. or len(p0) != 384**2:
        raise ValueError('Expected original 384² t0')
    xy0 = np.column_stack((p0['x'], p0['y'])).astype(float)
    # Match the production audit exactly, including sites on r0=0.5 or 2.
    r0 = np.hypot(xy0[:, 0]-3., xy0[:, 1]-3.)
    selected = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    foreign = r0 >= 2.75
    mass = p0['mass'][selected].astype(float)
    rho0 = p0['den'][selected].astype(float)
    volume0 = p0['volume'][selected].astype(float)
    initial_areas, _, _, safety0, poly0 = geometry(xy0, selected, foreign, True)
    initial_area_error = float(np.max(abs(initial_areas/volume0-1.)))
    if initial_area_error > 1.e-6:
        raise ValueError('Independent initial area reader mismatch')
    rows = []
    for checkpoint in checkpoints:
        p, meta = ordered(checkpoint)
        if not np.array_equal(p0['mass'], p['mass']) or not np.all(p['w2'] == 0):
            raise ValueError('Expected same fixed mass and W0')
        if not all(np.isfinite(p[k]).all() and np.all(p[k] > 0)
                   for k in ['mass', 'den', 'volume', 'pressure', 'ie']):
            raise ValueError('Nonfinite or nonpositive state')
        xy = np.column_stack((p['x'], p['y'])).astype(float)
        exact = orbit(xy0, meta['t'])
        restored = xy.copy()
        restored[foreign] = xy0[foreign]
        isolated = xy.copy()
        isolated[foreign] = exact[foreign]
        areas, support, edges, safety, polygons = geometry(xy, selected, foreign, True)
        restored_area, _, _, restored_safety, restored_polygons = geometry(restored, selected, foreign, True)
        isolated_area, _, _, _, _ = geometry(isolated, selected, foreign, False)
        exact_area, _, _, _, _ = geometry(exact, selected, foreign, False)
        native_rho = p['den'][selected].astype(float)
        native_volume = p['volume'][selected].astype(float)
        change = native_rho/rho0-1.
        pos = np.linalg.norm(xy[selected]-exact[selected], axis=1)
        entropy = (np.log(p['ie'][selected].astype(float)/p0['ie'][selected].astype(float))+
                   GAMMA_MINUS_ONE*np.log(native_volume/volume0))
        values = {'native_density_abs': abs(change), 'native_density_signed': change,
                  'independent_native_area_relative_abs': abs(areas/native_volume-1.),
                  'restored_exterior_density_abs': abs(mass/restored_area/rho0-1.),
                  'isolated_exterior_density_abs': abs(mass/isolated_area/rho0-1.),
                  'isolated_all_exact_density_abs': abs(mass/exact_area/rho0-1.),
                  'native_vs_isolated_all_exact_abs': abs(exact_area/native_volume-1.),
                  'entropy_log_abs': abs(entropy), 'entropy_log_signed': entropy,
                  'position_error': pos}
        bins = []
        for lo, hi in [(.5,.75),(.75,1.),(1.,2.),(.5,2.)]+[(x/20,(x+1)/20) for x in range(20,40)]:
            take = (r0[selected] >= lo) & (r0[selected] < hi)
            bins.append({'rlo':lo, 'rhi':hi, **row(take, mass, values)})
        contribution = mass*abs(change)
        support_map = {}
        for i,j in edges:
            support_map.setdefault(i, []).append(j)
        top = []
        for k in np.argsort(contribution)[-30:][::-1]:
            i = int(selected[k])
            top.append({'id':i, 'r0':float(r0[i]), 'r':float(np.linalg.norm(xy[i]-3.)),
                        'relative_density_change':float(change[k]),
                        'position_error':float(pos[k]), 'foreign_support_ids':support_map.get(i, [])})
        # Independent of Voronoi volumes: mass within fixed radial/spatial bins.
        radial_mass = []
        rt = np.hypot(xy[:, 0]-3., xy[:, 1]-3.)
        for lo, hi in [(1.,1.5),(1.5,1.75),(1.75,2.),(.5,2.)]:
            current = (rt >= lo) & (rt < hi)
            start = (r0 >= lo) & (r0 < hi)
            radial_mass.append({'rlo':lo, 'rhi':hi,
                'initial_mass':float(np.sum(p0['mass'][start],dtype=np.longdouble)),
                'current_mass':float(np.sum(p['mass'][current],dtype=np.longdouble)),
                'current_foreign_mass':float(np.sum(p['mass'][current & foreign],dtype=np.longdouble)),
                'foreign_generator_count':int(np.count_nonzero(current & foreign))})
        angle = np.arctan2(xy[selected,1]-3.,xy[selected,0]-3.)
        modes = {str(m):float(abs(np.sum(mass*change*np.exp(-1j*m*angle))/sum(mass)))
                 for m in range(1,9)}
        result = {'checkpoint':str(checkpoint), 'sha256':digest(checkpoint), 'time':meta['t'],
                  'disk_cells':len(selected), 'bins':bins,
                  'native_independent_area_max_relative':float(np.max(abs(areas/native_volume-1.))),
                  'actual_periodic_image_safety_min':safety,
                  'restored_periodic_image_safety_min':restored_safety,
                  'foreign_support_disk_cells':int(support.sum()), 'foreign_support_faces':len(edges),
                  'foreign_support_share_of_abs_density_norm':float(sum(contribution[support])/sum(contribution)),
                  'support_disk_density_abs':row(support,mass,values),
                  'no_foreign_support_disk_density_abs':row(~support,mass,values),
                  'foreign_support_initial_radius_min':float(min(r0[j] for _,j in edges)) if edges else None,
                  'largest_error_cells':top, 'fixed_annulus_mass':radial_mass,
                  'signed_density_fourier_amplitudes_m1_to8':modes}
        rows.append(result)
        print(json.dumps({'time':meta['t'],'disk':bins[3], 'foreign_support_cells':int(support.sum()),
                          'area_max_error':result['native_independent_area_max_relative']}),flush=True)
        if meta['t'] == 3.:
            render(output/'density_t3.png', poly0, polygons, restored_polygons,
                   rho0, native_rho, mass/restored_area)
        if result['native_independent_area_max_relative'] > 1.e-6:
            raise ValueError('Stored native disk volume disagrees with independent periodic geometry')
    report = {'initial':str(initial),'initial_sha256':digest(initial),
              'initial_area_max_relative':initial_area_error, 'initial_image_safety_min':safety0,
              'states':rows,'analysis_seconds':time.monotonic()-started,
              'scope':'Read-only saved states; actual/restored box configurations have certified complete periodic disk cells. Isolated exact-orbit counterfactuals intentionally use a different boundary/support configuration and are not alternate periodic solutions. No gas evolution, production modification or long-time proof.'}
    with (output/'support_comparison.json').open('x') as stream:
        json.dump(report,stream,indent=2,allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config',type=Path)
    args=parser.parse_args()
    config=json.loads(args.config.read_text())
    for path, checksum in config['input_hashes'].items():
        if digest(path)!=checksum:
            raise ValueError('Frozen input changed: '+path)
    analyze(Path(config['initial']),list(map(Path,config['checkpoints'])),Path.cwd())

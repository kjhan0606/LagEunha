"""Construct matched uniform/graded Kepler ICs; evolution is native CUDA.

The monitor below controls INITIAL mesh sampling only. It is not gas density,
Laguerre weight, viscosity, or a generator velocity during evolution.
"""
from pathlib import Path
import argparse
import hashlib
import json
import struct
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.spatial import Voronoi
from scipy.stats import qmc

from audit_snapshot import read_snapshot

L = 24.0
CENTER = 12.0
H0 = L / 256
RC = 0.25
EPS = 0.047
GAMMA = float(np.float32(1.666667))
PRESSURE = float(np.float32(1e-6))
SEED = 17012026
CVT_STEPS = 16
RADIAL_EDGES = [0, .047, .1, .2, .5, 1, 2, 4, 8, 18]


def monitor(r, graded):
    if not graded:
        return np.ones_like(r)
    s = np.clip((r - 2) / 2, 0, 1)
    taper = 1 - 10*s**3 + 15*s**4 - 6*s**5
    taper = np.where(r <= 2, 1, np.where(r >= 4, 0, taper))
    return 1 + 63*RC**2/(r*r + RC**2)*taper


def initial_points(graded):
    extra = quad(lambda r: 2*np.pi*r*(float(monitor(np.array(r), True))-1),
                 0, 4, points=[2], epsabs=1e-10, epsrel=1e-12)[0] / H0**2
    count = int(round(256**2 + (extra if graded else 0)))
    maximum = 64 if graded else 1
    sequence = qmc.Sobol(d=3, scramble=True, seed=SEED)
    selected = []
    total = 0
    generated = 0
    while total < count:
        block = sequence.random(65536)
        xy = block[:, :2] * L
        accept = block[:, 2]*maximum < monitor(np.linalg.norm(xy-CENTER, axis=1), graded)
        chosen = xy[accept]
        selected.append(chosen)
        total += len(chosen)
        generated += len(block)
    return np.vstack(selected)[:count].copy(), dict(
        seed=SEED, candidate_points=generated, count=count,
        target_integral=256**2+(extra if graded else 0),
        initialization='scrambled 3D Sobol rejection; same algorithm and seed')


def geometry(xy, statistics=False):
    """Full 3x3 images, plus a per-cell exclusion proof for farther images."""
    blocks = [xy] + [xy + np.array([dx, dy]) for dx in (-L, 0, L)
                    for dy in (-L, 0, L) if dx or dy]
    vor = Voronoi(np.vstack(blocks))
    volume = np.empty(len(xy))
    radius = np.empty(len(xy))
    apex, ends_a, ends_b, owner, central = [], [], [], [], []
    eigenratio = np.empty(len(xy)) if statistics else None
    width = np.empty(len(xy)) if statistics else None
    sides = np.empty(len(xy), dtype=int) if statistics else None
    for i in range(len(xy)):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError('unbounded periodic cell')
        v = vor.vertices[region] - CENTER
        local = v - (xy[i]-CENTER)
        cross = local[:, 0]*np.roll(local[:, 1], -1)-local[:, 1]*np.roll(local[:, 0], -1)
        if cross.sum() < 0:
            v, local = v[::-1], local[::-1]
            cross = local[:, 0]*np.roll(local[:, 1], -1)-local[:, 1]*np.roll(local[:, 0], -1)
        volume[i] = .5*cross.sum()
        radius[i] = np.max(np.linalg.norm(local, axis=1))
        if not (volume[i] > 0 and 2*radius[i] < L):
            raise ValueError('positive complete periodic cell proof failed')
        nxt = np.roll(v, -1, axis=0)
        edge = nxt-v
        signed = edge[:, 0]*(-v[:, 1])-edge[:, 1]*(-v[:, 0])
        roundoff = 128*np.finfo(float).eps*max(1, np.max(abs(v))**2)
        inside = bool(np.all(signed >= -roundoff))
        origin = np.zeros(2) if inside else xy[i]-CENTER
        for a, b in zip(v, nxt):
            da, db = a-origin, b-origin
            c = da[0]*db[1]-da[1]*db[0]
            if c < -roundoff:
                raise ValueError('nonconvex triangle fan')
            if c <= 0:
                continue
            apex.append(origin)
            ends_a.append(a)
            ends_b.append(b)
            owner.append(i)
            central.append(inside)
        if statistics:
            q = np.roll(local, -1, axis=0)
            centroid = np.sum(cross[:, None]*(local+q), axis=0)/(6*volume[i])
            xx = np.sum(cross*(local[:, 0]**2+local[:, 0]*q[:, 0]+q[:, 0]**2))/(12*volume[i])
            yy = np.sum(cross*(local[:, 1]**2+local[:, 1]*q[:, 1]+q[:, 1]**2))/(12*volume[i])
            xy2 = np.sum(cross*(2*local[:, 0]*local[:, 1]+local[:, 0]*q[:, 1]
                               +q[:, 0]*local[:, 1]+2*q[:, 0]*q[:, 1]))/(24*volume[i])
            cov = np.array([[xx, xy2], [xy2, yy]])-np.outer(centroid, centroid)
            ev = np.linalg.eigvalsh(cov)
            if ev[0] <= 0:
                raise ValueError('nonpositive polygon inertia')
            eigenratio[i] = ev[1]/ev[0]
            normal = np.column_stack((-edge[:, 1], edge[:, 0]))
            normal /= np.linalg.norm(normal, axis=1)[:, None]
            projection = local @ normal.T
            width[i] = np.min(np.ptp(projection, axis=0))
            sides[i] = len(v)
    if abs(volume.sum()-L*L) > 1e-8:
        raise ValueError('periodic total area failed')
    fan = (volume, radius, np.array(apex), np.array(ends_a), np.array(ends_b),
           np.array(owner), np.array(central))
    stats = None
    if statistics:
        pair = vor.ridge_points
        pair = pair[np.any(pair < len(xy), axis=1)] % len(xy)
        pair = np.sort(pair, axis=1)
        pair = np.unique(pair[pair[:, 0] != pair[:, 1]], axis=0)
        stats = dict(inertia_eigenvalue_ratio=eigenratio, minimum_width=width,
                     polygon_sides=sides, neighbor_pairs=pair)
    return fan, stats


def relax(xy, graded):
    history = []
    for iteration in range(CVT_STEPS):
        start = time.monotonic()
        fan, _ = geometry(xy)
        _, _, origin, a, b, owner, _ = fan
        da, db = a-origin, b-origin
        area = .5*(da[:, 0]*db[:, 1]-da[:, 1]*db[:, 0])
        mass = np.zeros(len(xy))
        moment = np.zeros((len(xy), 2))
        for aa, bb in ((1/6, 1/6), (2/3, 1/6), (1/6, 2/3)):
            point = origin+aa*da+bb*db
            # In 2D the CVT functional density is proportional to n^2.
            weight = area/3*monitor(np.linalg.norm(point, axis=1), graded)**2
            mass += np.bincount(owner, weights=weight, minlength=len(xy))
            for d in range(2):
                moment[:, d] += np.bincount(owner, weights=weight*point[:, d], minlength=len(xy))
        if not np.all(mass > 0):
            raise ValueError('nonpositive initialization CVT weight')
        updated = (moment/mass[:, None]+CENTER) % L
        shift = (updated-xy+L/2) % L-L/2
        normalized = np.linalg.norm(shift, axis=1)/np.sqrt(fan[0])
        xy = updated
        row = dict(iteration=iteration+1, max_shift_cell_length=float(normalized.max()),
                   rms_shift_cell_length=float(np.sqrt(np.mean(normalized**2))),
                   seconds=time.monotonic()-start)
        history.append(row)
        print('CVT', 'graded' if graded else 'uniform', row, flush=True)
    return xy, history


def density(r):
    shape = np.exp(-np.logaddexp(0., 8*(2-r))-np.logaddexp(0., 8*(r-8)))
    return .001+2*shape/(r+1e-6)


def integrate_density(n, fan, selected=None):
    volume, _, origin, a, b, owner, central = fan
    active = np.ones(len(owner), bool) if selected is None else selected[owner]
    indices = np.flatnonzero(active)
    result = np.zeros(len(volume))
    nodes, weights = leggauss(n)
    nodes, weights = (nodes+1)/2, weights/2
    u, v = nodes[None, :, None], nodes[None, None, :]
    wu, wv = weights[None, :, None], weights[None, None, :]
    logf0 = -np.logaddexp(0., 16.)-np.logaddexp(0., -64.)
    f0 = np.exp(logf0)
    for start in range(0, len(indices), 512):
        ix = indices[start:start+512]
        o, da, db = origin[ix], a[ix]-origin[ix], b[ix]-origin[ix]
        edge = (1-v[..., None])*da[:, None, None, :]+v[..., None]*db[:, None, None, :]
        points = o[:, None, None, :]+u[..., None]*edge
        radius = np.linalg.norm(points, axis=-1)
        # Subtract the known central radial term on EVERY cell, not only
        # the polygon containing the origin. Its exact area integral is a
        # boundary integral by div(x G(r))=1/(r+a),
        # G(r)=[r-a log(1+r/a)]/r^2. This retains the original density.
        logf = -np.logaddexp(0., 8*(2-radius))-np.logaddexp(0., 8*(radius-8))
        rho = .001+2*f0*np.expm1(logf-logf0)/(radius+1e-6)
        value = np.sum(wu*wv*u*rho, axis=(1, 2))
        cross = da[:, 0]*db[:, 1]-da[:, 1]*db[:, 0]
        result += np.bincount(owner[ix], weights=cross*value, minlength=len(volume))
        line = ((1-nodes)[None, :, None]*a[ix, None, :]
                +nodes[None, :, None]*b[ix, None, :])
        rr = np.linalg.norm(line, axis=-1)
        z = rr/1e-6
        primitive = z-np.log1p(z)
        small = z < .01
        if small.any():
            zz = z[small]
            # Alternating convergent series avoids cancellation at r=0.
            primitive[small] = zz**2*(.5+zz*(-1/3+zz*(.25+zz*(-.2+zz*(1/6+zz*(-1/7+zz/8))))))
        green = primitive*1e-6/(rr*rr)
        boundary = 2*f0*(a[ix, 0]*b[ix, 1]-a[ix, 1]*b[ix, 0])*np.sum(weights[None, :]*green, axis=1)
        result += np.bincount(owner[ix], weights=boundary, minlength=len(volume))
    return result


def resolve_near_origin_edges(fan):
    """Resolve the fixed analytic radial boundary primitive, not the mesh."""
    volume, radius, origin, a, b, owner, central = fan
    oo, aa, bb, ids, cc = [], [], [], [], []
    for o, x, y, i, c in zip(origin, a, b, owner, central):
        edge = y-x
        length2 = float(edge @ edge)
        t = float(np.clip(-(x @ edge)/length2, 0, 1))
        distance = float(np.linalg.norm(x+t*edge))
        step = (distance+1e-6)/np.sqrt(length2)
        cuts = [0., 1.]
        if step < .25:
            cuts.append(t)
            for k in range(32):
                width = step*2**k
                if width >= 2:
                    break
                for point in (t-width, t+width):
                    if 0 < point < 1:
                        cuts.append(point)
        cuts = sorted(set(cuts))
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            oo.append(o)
            aa.append(x+lo*edge)
            bb.append(x+hi*edge)
            ids.append(i)
            cc.append(c)
    return volume, radius, np.array(oo), np.array(aa), np.array(bb), np.array(ids), np.array(cc)


def radial_statistics(xy, volume, stats, graded):
    r = np.linalg.norm(xy-CENTER, axis=1)
    h = np.sqrt(volume)
    pairs = stats['neighbor_pairs']
    adjacent_ratio = np.maximum(h[pairs[:, 0]], h[pairs[:, 1]])/np.minimum(h[pairs[:, 0]], h[pairs[:, 1]])
    rows = []
    for lo, hi in zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:]):
        selected = (r >= lo) & (r < hi)
        if not selected.any():
            continue
        neighbors = selected[pairs[:, 0]] | selected[pairs[:, 1]]
        rows.append(dict(Rlo=lo, Rhi=hi, n=int(selected.sum()),
                         median_h=float(np.median(h[selected])),
                         median_h_over_target=float(np.median(h[selected]*np.sqrt(monitor(r[selected], graded))/H0)),
                         neighbor_h_ratio_max=float(adjacent_ratio[neighbors].max()),
                         inertia_eigenvalue_ratio_p99=float(np.quantile(stats['inertia_eigenvalue_ratio'][selected], .99)),
                         inertia_eigenvalue_ratio_max=float(stats['inertia_eigenvalue_ratio'][selected].max()),
                         minimum_width=float(stats['minimum_width'][selected].min()),
                         neighbor_count_min=int(stats['polygon_sides'][selected].min()),
                         neighbor_count_max=int(stats['polygon_sides'][selected].max())))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--historical-root', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    layout, _ = read_snapshot(args.historical_root/'core_refinement/coarse_integrated/khout.000001.dat')
    params = (args.historical_root/'core_refinement/coarse_integrated/params.dat').read_text()
    manifest = dict(purpose='user-approved matched smooth radial-resolution comparison; not long-run acceptance',
                    h_outer=H0, r_core=RC, eps=EPS, gamma=GAMMA, pressure=PRESSURE,
                    monitor='1+63*rc^2/(r^2+rc^2)*C2_taper(r;2,4)',
                    cvt_steps=CVT_STEPS, initialization_only=True, dynamic_steering=0,
                    weights=0, Re='infinity', individual_masses_fixed_during_evolution=True,
                    initial_velocity='same original pointwise Plummer circular velocity at each generator',
                    output_policy='two ICs and compact geometry arrays, approximately 90 MiB; no evolution dumps here',
                    cases={})
    for name, graded in (('uniform', False), ('graded', True)):
        start = time.monotonic()
        xy, sampling = initial_points(graded)
        xy, history = relax(xy, graded)
        dest = args.output/name
        dest.mkdir()
        np.savez_compressed(dest/'initialized_points.npz', xy=xy)
        (dest/'initialization.json').write_text(json.dumps(dict(sampling=sampling, cvt_history=history), indent=2)+'\n')
        fan, stats = geometry(xy, statistics=True)
        volume = fan[0]
        fan = resolve_near_origin_edges(fan)
        old, mass = integrate_density(4, fan), integrate_density(8, fan)
        error = abs(mass-old)/mass
        refinements = []
        for order in (16, 32, 64):
            bad = error > 1e-8
            if not bad.any():
                break
            revised = integrate_density(order, fan, bad)
            error[bad] = abs(revised[bad]-mass[bad])/revised[bad]
            mass[bad] = revised[bad]
            refinements.append(dict(order=order, cells=int(bad.sum())))
        if not np.all(mass > 0) or error.max() > 1e-8:
            np.savez_compressed(dest/'failed_integration.npz', mass=mass, error=error)
            print('INTEGRATION_FAILED', name, 'max_relative_change', error.max(),
                  'worst_cell', int(np.argmax(error)), 'radius', np.linalg.norm(xy[np.argmax(error)]-CENTER), flush=True)
            raise ValueError('original density integration did not converge')
        raw = np.zeros(len(xy)*layout.dtype.itemsize, dtype=np.uint8)
        p = np.ndarray(len(xy), dtype=layout.dtype, buffer=raw)
        radius = np.linalg.norm(xy-CENTER, axis=1)
        omega = (radius*radius+EPS*EPS)**(-.75)
        p['u4if'] = np.arange(len(xy), dtype=np.uint64)
        p['x'], p['y'] = xy[:, 0], xy[:, 1]
        p['mass'], p['volume'], p['den'] = mass, volume, mass/volume
        p['pressure'], p['ie'] = PRESSURE, PRESSURE*volume/(GAMMA-1)
        p['vx'], p['vy'] = -(xy[:, 1]-CENTER)*omega, (xy[:, 0]-CENTER)*omega
        path = dest/'khout.000001.dat'
        with path.open('wb') as f:
            f.write(struct.pack('<4i', -1, len(p), layout.dtype.itemsize, 1))
            f.write(raw.tobytes())
            f.write(struct.pack('<ddfiiddff', 0., 1e-7, 1e-7, 256, 256, L, L, 0., 0.))
        (dest/'params.dat').write_text(params)
        np.savez_compressed(dest/'geometry.npz', xy=xy, radius=radius, volume=volume,
                            mass=mass, quad_error=error, **stats)
        kinetic = .5*mass*(p['vx']**2+p['vy']**2)
        potential = -mass/np.sqrt(radius**2+EPS**2)
        angular = mass*((xy[:, 0]-CENTER)*p['vy']-(xy[:, 1]-CENTER)*p['vx'])
        row = dict(sampling=sampling, cvt_history=history, n=len(p),
                   total_mass=float(mass.sum()), total_area=float(volume.sum()),
                   total_internal_energy=float(p['ie'].astype(float).sum()),
                   total_kinetic_energy=float(kinetic.sum()), total_potential_energy=float(potential.sum()),
                   total_angular_momentum=float(angular.sum()),
                   min_generator_radius=float(radius.min()),
                   max_relative_mass_quadrature_change=float(error.max()),
                   quadrature_refinements=refinements,
                   radial_mesh=radial_statistics(xy, volume, stats, graded),
                   IC_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                   preparation_seconds=time.monotonic()-start)
        manifest['cases'][name] = row
        (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
        print('IC_COMPLETE', name, 'N', len(p), 'mass', mass.sum(), flush=True)
    a, b = manifest['cases'].values()
    difference = abs(a['total_mass']-b['total_mass'])/a['total_mass']
    manifest['uniform_graded_mass_relative_difference'] = difference
    if difference > 1e-8:
        raise ValueError('matched continuous mass contract failed')
    manifest['IC_preparation_pass'] = True
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    main()

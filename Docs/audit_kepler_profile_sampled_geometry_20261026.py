"""Independent Voronoi geometry gate for profile-sampled Kepler glass."""

import argparse
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import Voronoi

EPS = .01175


def glass(path):
    raw = path.read_bytes()
    n, lx, ly = struct.unpack_from('<idd', raw)
    if not (0 < n <= 384*384 and (lx, ly) == (6., 6.) and len(raw) == 20+16*n):
        raise ValueError('Invalid glass file')
    x = np.frombuffer(raw, dtype='<f8', count=n, offset=20)
    y = np.frombuffer(raw, dtype='<f8', count=n, offset=20+8*n)
    return np.column_stack((x, y))


def areas(points, selected):
    vor = Voronoi(points, qhull_options='Qbb Qc Qz')
    result = np.empty(len(selected))
    for k, i in enumerate(selected):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f'Unbounded selected cell {i}')
        v = vor.vertices[region]
        result[k] = .5*abs(np.dot(v[:, 0], np.roll(v[:, 1], -1))-
                              np.dot(v[:, 1], np.roll(v[:, 0], -1)))
    if not np.all(np.isfinite(result) & (result > 0)):
        raise ValueError('Invalid Voronoi area')
    return result


def orbit(points, t):
    dx, dy = points[:, 0]-3., points[:, 1]-3.
    r2 = dx*dx+dy*dy
    phase = t*(r2+EPS*EPS)**(-.75)
    return np.column_stack((3.+dx*np.cos(phase)-dy*np.sin(phase),
                            3.+dx*np.sin(phase)+dy*np.cos(phase)))


def weighted_mean(x, w):
    return float(np.sum(x*w, dtype=np.longdouble)/np.sum(w, dtype=np.longdouble))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=Path)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    ref, cand = glass(args.reference), glass(args.candidate)
    rr = np.hypot(ref[:, 0]-3., ref[:, 1]-3.)
    rc = np.hypot(cand[:, 0]-3., cand[:, 1]-3.)
    ir = np.flatnonzero((rr>=.5) & (rr<2.))
    ic = np.flatnonzero((rc>=.5) & (rc<2.))
    same_disk = bool(np.array_equal(ref[ir], cand[ic]))
    if not same_disk:
        raise ValueError('Candidate did not preserve q=2 disk generators')
    ar0 = areas(ref, ir)
    ac0 = areas(cand, ic)
    mass = 1.01*ac0
    area_difference = abs(ac0/ar0-1)
    radial = rc[ic]
    bins = []
    for lo, hi in ((.5,.6),(.6,1.),(1.,1.9),(1.9,2.),(.5,2.)):
        use = (radial>=lo)&(radial<hi)
        err = area_difference[use]
        bins.append({'rlo':lo,'rhi':hi,'cells':int(np.sum(use)),
                     'initial_area_relative_abs_mass_weighted':weighted_mean(err,mass[use]),
                     'initial_area_relative_p95':float(np.quantile(err,.95)),
                     'initial_area_relative_max':float(np.max(err)),
                     'candidate_cell_mass_cv':float(np.std(mass[use])/np.mean(mass[use]))})
    times=[]
    for t in (.20626460064652694, 2.221441469079183, 600.):
        act = areas(orbit(cand,t), ic)
        change = abs(ac0/act-1)
        rows=[]
        for lo,hi in ((.5,.75),(.75,1.),(1.,2.),(.5,2.)):
            use=(radial>=lo)&(radial<hi)
            rows.append({'rlo':lo,'rhi':hi,'density_change_mass_weighted':weighted_mean(change[use],mass[use]),
                         'p95':float(np.quantile(change[use],.95))})
        times.append({'time':t,'bins':rows})
    result={'reference':str(args.reference),'candidate':str(args.candidate),
            'reference_sites':len(ref),'candidate_sites':len(cand),
            'disk_generators_bitwise_preserved':same_disk,
            'disk_initial_area_comparison':bins,'ballistic_candidate':times,
            'scope':'Exact-orbit Voronoi geometry only; not a native hydro result'}
    with args.output.open('x') as f:
        json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()

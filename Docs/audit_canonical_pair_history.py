"""Read-only pair census of existing canonical checkpoints; run on lageunha.

Sparse snapshots cannot prove continuous sticking or identify individual forces.
The reference is the exact circular orbit in the configured external potential,
not a replacement simulation or a modification of the saved material states.
"""
import argparse
import hashlib
import json
import struct
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def gravity_coefficient(r, rc=.5):
    u = r / rc
    enclosed = np.ones_like(r)
    lo, mid = u <= .5, (u > .5) & (u < 1)
    q = u[lo]
    enclosed[lo] = q**3 * (32/3 - 192/5*q*q + 32*q**3)
    q = u[mid]
    enclosed[mid] = -1/15 + q**3*(64/3-48*q+192/5*q*q-32/3*q**3)
    return enclosed / np.maximum(r, 1e-30)**3


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    base = Path('/gpfs/kjhan/LagForce')
    layout_path = base/'Codex_DPQ_backend_parity_20261017/layout.json'
    layout = json.loads(layout_path.read_text())
    names = ['x', 'y', 'vx', 'vy', 'mass', 'volume', 'den', 'pressure', 'ie']
    dtype = np.dtype({'names': ['id']+names,
                      'formats': ['<u8']+['<f8']*len(names),
                      'offsets': [layout['u4if']]+[layout[k] for k in names],
                      'itemsize': layout['size']})

    def read(path):
        raw = path.read_bytes()
        magic, n, size, stress = struct.unpack_from('<iiii', raw)
        assert magic == -1 and stress == 1 and size == dtype.itemsize
        assert len(raw) == 16+n*size+52
        a = np.frombuffer(raw, dtype=dtype, count=n, offset=16)
        ids = a['id'] & layout['id_mask']
        order = np.argsort(ids)
        assert np.array_equal(ids[order], np.arange(n))
        a = a[order]
        assert all(np.isfinite(a[k]).all() for k in names)
        assert all(np.all(a[k] > 0) for k in ['mass','volume','den','pressure','ie'])
        return a, struct.unpack_from('<d',raw,16+n*size)[0], hashlib.sha256(raw).hexdigest()

    initial_path = base/'Codex_canonical_initial_20261017/run/khout.000000.dat'
    initial, t0, initial_sha = read(initial_path)
    assert t0 == 0
    x0 = np.column_stack((initial['x']-3, initial['y']-3))
    v0 = np.column_stack((initial['vx'], initial['vy']))
    r0 = np.linalg.norm(x0, axis=1)
    omega = np.sqrt(gravity_coefficient(r0))
    vref0 = omega[:,None]*np.column_stack((-x0[:,1], x0[:,0]))
    active = r0 < 2.4  # away from the frozen boundary at 2.5
    initial_velocity_error = float(np.max(np.linalg.norm((vref0-v0)[active],axis=1)))
    assert initial_velocity_error < 1e-12
    regions = {
        'inner_background': (r0 < .5) & (initial['den'] < .5),
        'disk_inner_edge': (r0 < .6) & (initial['den'] > .5),
        'disk_core': (r0 >= .6) & (r0 < 1.9) & (initial['den'] > .5),
    }
    run = base/'Codex_HG_constrained_20261017/omp8_t20'
    paths = [initial_path]+[run/f'khout.{i:06d}.dat' for i in [184,368,656,1267,2305,4114,6248]]
    states, positions, velocities, reference, refvel, times, sources = [], [], [], [], [], [], []
    candidates = set()
    census = []
    h0 = np.sqrt(initial['volume'])
    for path in paths:
        a, t, sha = read(path)
        assert np.array_equal(a['mass'], initial['mass'])
        x = np.column_stack((a['x']-3,a['y']-3))
        v = np.column_stack((a['vx'],a['vy']))
        c,s = np.cos(omega*t),np.sin(omega*t)
        xr = np.column_stack((c*x0[:,0]-s*x0[:,1],s*x0[:,0]+c*x0[:,1]))
        vr = omega[:,None]*np.column_stack((-xr[:,1],xr[:,0]))
        dd,jj = cKDTree(x).query(x,k=2,workers=1)
        dr,jr = cKDTree(xr).query(xr,k=2,workers=1)
        near, near_ref = dd[:,1], dr[:,1]
        j = jj[:,1]
        scale = np.sqrt(.5*(initial['volume']+initial['volume'][j]))
        # Keep each inner particle's nearest pair at every saved time,
        # not only the pairs surviving to the last checkpoint.
        for i in np.flatnonzero(r0 < .7):
            candidates.add(tuple(sorted((int(i),int(j[i])))))
        row = {'time': t, 'regions': {}}
        for name, mask in regions.items():
            q = near[mask]/scale[mask]
            qr = near_ref[mask]/np.sqrt(.5*(initial['volume'][mask]+initial['volume'][jr[mask,1]]))
            err = a['pressure']/initial['pressure']-1
            row['regions'][name] = {
                'cells': int(mask.sum()),
                'actual_gap_over_initial_h_quantiles': np.quantile(q,[0,.01,.1,.5]).tolist(),
                'circular_gap_over_initial_h_quantiles': np.quantile(qr,[0,.01,.1,.5]).tolist(),
                'actual_fraction_below_0p2': float(np.mean(q < .2)),
                'circular_fraction_below_0p2': float(np.mean(qr < .2)),
                'pressure_error_mass_RMS': float(np.sqrt(np.average(err[mask]**2,weights=a['mass'][mask]))),
                'radial_error_mass_RMS': float(np.sqrt(np.average((np.linalg.norm(x[mask],axis=1)-r0[mask])**2,weights=a['mass'][mask]))),
            }
        census.append(row)
        states.append(a); positions.append(x); velocities.append(v)
        reference.append(xr); refvel.append(vr); times.append(t)
        sources.append({'path':str(path),'time':t,'sha256':sha})
    pairs = np.asarray(sorted(candidates), dtype=int)
    i,j = pairs.T
    pair_h = np.sqrt(.5*(initial['volume'][i]+initial['volume'][j]))
    def histories(xx,vv):
        dx = np.asarray(xx)[:,j]-np.asarray(xx)[:,i]
        dv = np.asarray(vv)[:,j]-np.asarray(vv)[:,i]
        d = np.linalg.norm(dx,axis=2)
        assert np.all(d > 0)
        radial = np.sum(dx*dv,axis=2)/d
        tangential = (dx[:,:,0]*dv[:,:,1]-dx[:,:,1]*dv[:,:,0])/d
        return d, radial, tangential
    d,vd,vt = histories(positions,velocities)
    dc,vdc,vtc = histories(reference,refvel)
    p = np.asarray([a['pressure'][pairs]/initial['pressure'][pairs] for a in states])
    volume = np.asarray([a['volume'][pairs]/initial['volume'][pairs] for a in states])
    q = d/pair_h
    close = q < .2
    # Count consecutive sampled states; explicitly not continuous residence.
    streak = np.zeros(len(pairs),dtype=int); longest = streak.copy()
    for row in close:
        streak = np.where(row,streak+1,0); longest = np.maximum(longest,streak)
    ranked = sorted(range(len(pairs)),key=lambda k:(-int(longest[k]),float(q[:,k].min())))[:16]
    selected = []
    for k in ranked:
        selected.append({'ids':pairs[k].tolist(),'initial_radii':r0[pairs[k]].tolist(),
                         'max_consecutive_saved_close_states':int(longest[k]),
                         'gap_over_initial_h':q[:,k].tolist(),
                         'same_pair_circular_gap_over_initial_h':(dc[:,k]/pair_h[k]).tolist(),
                         'separation_rate':vd[:,k].tolist(),
                         'transverse_relative_speed':vt[:,k].tolist(),
                         'pressure_over_initial':p[:,k].tolist(),
                         'volume_over_initial':volume[:,k].tolist()})
    np.savez_compressed(args.output/'pair_history.npz',pairs=pairs,times=times,h0=pair_h,
                        separation=d,circular_separation=dc,separation_rate=vd,
                        circular_separation_rate=vdc,transverse_speed=vt,
                        circular_transverse_speed=vtc,pressure_ratio=p,volume_ratio=volume)
    result = {'scope':'Existing original fixed-HG canonical checkpoints, no new evolution or model changes.',
              'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'layout':str(layout_path),'initial_sha256':initial_sha,
              'initial_circular_velocity_max_error':initial_velocity_error,
              'normalization':'h_ij=sqrt((V_i(t0)+V_j(t0))/2); fixed material scale.',
              'reference':'Exact circular point orbits in spline external potential rc=.5, GM=1, center=(3,3).',
              'sources':sources,'census':census,'candidate_pairs':len(pairs),
              'pairs_close_at_two_or_more_consecutive_saved_states':int(np.sum(longest>=2)),
              'selected_pairs':selected,
              'limitations':['Snapshots spaced about 2 time units; repeated small gaps do not prove continuous sticking.',
                             'Nearest geometric neighbor need not share a Voronoi face.',
                             'No pressure/TCF/HG force decomposition yet; correlation is not causation.',
                             'Inner candidate selection r(t0)<.7; census also covers disk core.',
                             'Euclidean distances used for central regions, away from periodic boundaries.']}
    (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'output':str(args.output),'candidate_pairs':len(pairs),
                      'repeated_close_pairs':result['pairs_close_at_two_or_more_consecutive_saved_states'],
                      'census':census,'top_pair':selected[0]},indent=2))


if __name__ == '__main__':
    main()

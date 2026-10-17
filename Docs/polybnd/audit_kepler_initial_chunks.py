"""Resume read-only reconstruction audits of completed Kepler IC chunks.

Missing integration chunks end the audited prefix; they are not failures.
This validates a specified density quadrature, not an evolved solution.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import numpy as np
import moment_ale_evolve as model
from moment_ale_evolve import decode
from factor_moment_evolve import rebuild
from kepler_disk_reference import polar_rule, ring_sites


def digest(data):
    return hashlib.sha256(data).hexdigest()


def atomic(path, data):
    p = Path(path)
    tmp = Path(str(p)+'.tmp')
    tmp.write_text(json.dumps(data, indent=2)+'\n')
    tmp.replace(p)


def run(root, order=16, max_chunks=None):
    root = Path(root)
    if max_chunks is not None and max_chunks<1:
        raise ValueError('positive chunk limit required')
    from density_quadrature import validate_order
    order = validate_order(order)
    os.environ['LAGFORCE_DENSITY_QUADRATURE_ORDER'] = str(order)
    identity_data = (root/'identity.json').read_bytes()
    identity = json.loads(identity_data)
    geometry_data = (root/'geometry.json').read_bytes()
    cells = json.loads(geometry_data)
    sites = ring_sites(identity['n'], require_count=identity['n']==192)
    if digest(sites.tobytes()) != identity['sites_sha256'] or len(cells)!=len(sites):
        raise ValueError('initial geometry identity mismatch')
    audit_identity = dict(initial_identity_sha256=digest(identity_data),
                          geometry_sha256=digest(geometry_data), density_order=order,
                          source_manifest_sha256=digest(Path(model.__file__).with_name('SHA256SUMS').read_bytes()),
                          script_sha256=digest(Path(__file__).read_bytes()))
    summaries = []
    for count, first in enumerate(range(0,len(cells),identity['chunk_size'])):
        if max_chunks is not None and count>=max_chunks:
            break
        path = root/f'cells_{first:06d}.json'
        if not path.exists():
            break
        data = path.read_bytes()
        chunk = json.loads(data)
        last = min(first+identity['chunk_size'],len(cells))
        geometry = digest(json.dumps([sites[first:last].tolist(),cells[first:last]]).encode())
        if chunk['first']!=first or chunk['geometry_sha256']!=geometry or len(chunk['rows'])!=last-first:
            raise ValueError('initial chunk geometry mismatch')
        key = dict(identity=audit_identity, chunk_sha256=digest(data))
        cache = Path(f'chunk_{first:06d}.json')
        if cache.exists():
            saved = json.loads(cache.read_text())
            if saved['key']!=key:
                raise ValueError('reconstruction audit cache mismatch')
        else:
            rows = []
            for i, item in enumerate(chunk['rows'], first):
                m = item['mass']
                state = rebuild(np.r_[item['row'],np.zeros(4)],np.array([m]),np.array([item['factor']]))
                raw = np.r_[m,state[3:20]]
                decoded = decode(np.asarray(cells[i]),raw)
                den = decoded['density']
                integrals = []
                for q in (16,24):
                    x,w = polar_rule(np.asarray(cells[i]),q,radii=(.25,.5,2.,2.5))
                    integrals.append(float(w@den.value(x)))
                rows.append(dict(cell=i, pressure=float(decoded['P']),
                    density_fit_residual=den.residual,
                    independent_mass_relative=abs(integrals[1]-m)/m,
                    independent_mass_quadrature_difference=abs(integrals[1]-integrals[0])/m))
            saved = dict(key=key,rows=rows)
            atomic(cache,saved)
        summaries.extend(saved['rows'])
        print('audited initial cells',len(summaries),'/',len(cells),flush=True)
    report = dict(audited_cells=len(summaries),total_cells=len(cells),density_order=order,
                  covers_all_initial_cells=len(summaries)==len(cells),full_validation=False,
                  max_independent_mass_relative=max((r['independent_mass_relative'] for r in summaries),default=None),
                  max_quadrature_difference=max((r['independent_mass_quadrature_difference'] for r in summaries),default=None),
                  min_pressure=min((r['pressure'] for r in summaries),default=None),
                  scope='Initial per-cell reconstruction only; no topology-wide closure, potential enrichment, or evolution validation.')
    atomic('report.json',report)
    print(json.dumps(report),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('root')
    parser.add_argument('--order',type=int,default=16)
    parser.add_argument('--max-chunks',type=int)
    args = parser.parse_args()
    run(args.root,args.order,args.max_chunks)

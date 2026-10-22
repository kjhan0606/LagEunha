"""Common Eulerian sampling prevents different cell-center bin coverage bias."""
from pathlib import Path
import json
import numpy as np
from scipy.spatial import cKDTree
from audit_snapshot import read_snapshot
from prepare_radial_resolution_v1 import density

root = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
base = root/'radial_resolution_20261021_v2'
result = dict(scope='identical fixed Eulerian pixel samples of W0 piecewise-constant density vs original stationary density',
              formula='sum(abs(rho_num-rho_exact))/sum(rho_exact) on common spatial bin',
              caveat='pixel integration estimate; not a rigorous continuum norm bound', cases={})
grids = []
for extent, n, edges in ((.5, 1024, [0,.047,.1,.2,.5]), (9., 1024, [.5,1,2,4,8])):
    line = (np.arange(n)+.5)/n*2*extent-extent
    x, y = np.meshgrid(line, line)
    radius = np.hypot(x, y).ravel()
    query = np.column_stack((x.ravel()+12, y.ravel()+12))
    exact = density(radius)
    grids.append((query, radius, exact, edges))
for label in ('uniform', 'graded'):
    native = json.loads((base/'native_v2'/label/'manifest.json').read_text())
    final = next(x['checkpoint'] for x in native['stages'] if x['label']=='t0.05' and x['exit_code']==0)
    row = {}
    for tag, path in (('initial', base/label/'khout.000001.dat'), ('t005', Path(final))):
        p, meta = read_snapshot(path)
        tree = cKDTree(np.column_stack((p['x'], p['y'])) % 24, boxsize=24)
        rho = p['mass']/p['volume']
        rows = []
        for query, radius, exact, edges in grids:
            indices = tree.query(query, workers=1)[1]
            error = abs(rho[indices]-exact)
            for lo, hi in zip(edges[:-1], edges[1:]):
                sel = (radius >= lo) & (radius < hi)
                rows.append(dict(Rlo=lo,Rhi=hi,pixels=int(sel.sum()),
                    rho_L1_relative=float(error[sel].sum()/exact[sel].sum())))
        row[tag] = dict(time=meta['t'], radial=rows)
    result['cases'][label] = row
(base/'eulerian_density_t005.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))

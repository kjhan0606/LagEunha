"""Compare preserved moving-cell endpoints; empirical temporal order only."""
import hashlib
import json
from pathlib import Path
import numpy as np

root=Path('/gpfs/kjhan/LagForce/Codex_coupled_relaxation_evolution_20261018')
cases=[]
for alpha in (0,1):
    states=[];metadata=[];sources=[]
    for dt,folder in ((.001,'results'),(.0005,'results'),(.00025,'results_refinement')):
        path=root/folder/f'alpha{alpha}_dt{dt:g}.npz'
        data=np.load(path)
        meta=json.loads(path.with_suffix('.json').read_text())
        assert meta['status']=='completed' and float(data['time'])==.01
        if states:
            assert np.array_equal(initial,data['initial']) and np.array_equal(mass,data['mass'])
        initial=data['initial'].copy();mass=data['mass'].copy()
        states.append(data['state'].copy());metadata.append(meta)
        sources.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    metrics={}
    selectors={'all_state':slice(None),'cell_moments':None}
    for name,selector in selectors.items():
        values=states if selector is not None else [s[:-4].reshape(-1,20)[:,3:] for s in states]
        diffs=[float(np.max(np.abs(values[i]-values[i+1]))) for i in (0,1)]
        ratio=diffs[0]/diffs[1]
        metrics[name]=dict(max_differences=diffs,ratio=ratio,observed_order=float(np.log2(ratio)))
    cases.append(dict(alpha=alpha,sources=sources,metrics=metrics,
        finest=metadata[-1],same_topology=all(m['final']['topology']==metadata[0]['final']['topology'] for m in metadata)))
print(json.dumps(dict(scope=__doc__,cases=cases,
    limitation='Three time steps on a hot 4x4 shear patch do not validate residual flux closure, spatial accuracy, or long Kepler evolution.'),indent=2))

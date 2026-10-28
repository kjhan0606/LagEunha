"""Spawn-safe rank-owned preparation using the unchanged checked cell kernels.

Only ownership/transport changes. All stencils, bounds and fit parameters are
global. A worker never loads MPI or initializes CUDA, nor forks another worker.
"""
import json, os, resource
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from full_moment_state_io import load as load_state

_key=None
_density=None
_velocity=None

def initialize():
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    assert os.environ['OPENBLAS_NUM_THREADS']=='1'
    assert os.environ['OMP_NUM_THREADS']=='1'

def load(descriptor):
    global _key,_density,_velocity
    descriptor=str(descriptor)
    if _key==descriptor:return _density,_velocity
    import full_disk_density_memory as d
    import full_disk_velocity_memory as v
    with Path(descriptor).open() as f:ctx=json.load(f)
    initial=load_state(ctx['state'])
    root=Path(ctx['geometry'])
    with np.load(root/'mesh.npz') as f:mesh={k:f[k] for k in f.files}
    native=json.loads((root/'full_exact_supports.json').read_text())
    mass,ids,Psi=(initial[k] for k in ('mass','ids','Psi'))
    n=len(ids);a=initial['state'][:-4].reshape(n,20)
    assert n==ctx['ncells'] and np.array_equal(mesh['ids'],ids)
    assert np.all(np.diff(ids)>0), 'Registry must be stable-ID sorted'
    assert np.array_equal(mesh['sites'],a[:,:2])
    assert np.array_equal(mesh['weights'],a[:,2])
    assert not initial['frozen'].any() and initial['psi_constrained_mask'].all()
    assert Psi.dtype==np.dtype(np.longdouble)
    raw=np.column_stack([mass,a[:,3:]])
    areas=mesh['moments'][:,0];centers=mesh['moments'][:,1:3]/areas[:,None]
    lookup={int(g):i for i,g in enumerate(ids)};neighbors=[set() for _ in ids]
    for x,y in mesh['face_keys']:
        i,j=lookup[int(x)],lookup[int(y)];neighbors[i].add(j);neighbors[j].add(i)
    bounds=json.loads(Path(ctx['bounds']).read_text()) if ctx.get('bounds') else None
    if bounds is not None:
        assert len(bounds)==n and all(row['status']=='passed' for row in bounds)
        bounds=[row['bounds'] for row in bounds]
    density=SimpleNamespace(mesh=mesh,cells=native['cells'],R=float(mesh['radius']),
        areas=areas,mass=mass,ids=ids,raw=raw,centers=centers,means=mass/areas,
        bounds=bounds,neighbors=neighbors,Psi=Psi,independent_q48=ctx['q48'])
    velocity=None
    if ctx.get('fits'):
        rows=json.loads(Path(ctx['fits']).read_text())
        assert [r['id'] for r in rows]==ids.tolist()
        assert all(r['status']=='passed' for r in rows)
        v.model.GAMMA=float(np.float32(1.666667))
        states=[v.model.decode(None,r,density=SimpleNamespace(mass=float(r[0]),
            area=float(row['parameter']['area']),origin=np.array(row['parameter']['origin'])))
            for r,row in zip(raw,rows)]
        factors=initial['factors']
        for state,L in zip(states,factors):state['R']=L@L.T
        centers=np.array([s['c'] for s in states]);means=np.array([s['v'] for s in states])
        offsets=np.r_[0,np.cumsum([r['nodes24'] for r in rows])].astype(np.int64)
        # Full registry offsets are shifted, not repartitioned or reordered.
        offsets=offsets-offsets[ctx['lo']]
        velocity=SimpleNamespace(fit_rows=rows,raw=raw,cells=native['cells'],
            R=float(mesh['radius']),mass=mass,factors=factors,centers=centers,
            means=means,second=raw[:,3:7].reshape(n,2,2)/mass[:,None,None],
            tree=v.cKDTree(centers),offsets=offsets,ids=ids,positions=a[:,:2],
            states=states,stage_key=ctx['stage_key'])
        ng=int(offsets[ctx['hi']])
        v.nodes=np.memmap(ctx['gravity'],dtype=np.float64,mode='r+',shape=(ng,7))
        v.cache=v.VolumeCache(1);v.cache.begin_stage(ctx['stage_key'])
    d._context=density;v._context=velocity
    _key=descriptor;_density=d;_velocity=v
    return d,v

def task(arg):
    descriptor,command,i=arg
    d,v=load(descriptor)
    if command==2:return d.bounds_worker(i)
    if command==3:return d.worker(i)
    if command==4:return v.worker(i)
    raise ValueError('Unknown cell preparation command')

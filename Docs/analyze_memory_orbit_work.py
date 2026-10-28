"""Recompute completed diagnostic branches without launching simulations.

Run on lageunha with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1.
Print JSON; source outputs are read-only. HG and RK remainder stay combined.
"""
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read, sha
from audit_canonical_pair_history import gravity_coefficient

root = Path('/gpfs/kjhan/LagForce/Codex_memory_orbit_work_control_20261018')
manifest = json.loads((root/'manifest.json').read_text())
assert sha(Path(manifest['binary'])) == manifest['sha256']
results = {}

def departure(a):
    x = np.column_stack((a['x']-3, a['y']-3))
    r = np.linalg.norm(x, axis=1)
    vc = np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack((-x[:,1], x[:,0]))
    u = np.column_stack((a['vx'], a['vy']))-vc
    active = r < 2.5
    return a['id'][active], float((.5*a['mass'][active]*np.sum(u[active]**2, axis=1)).sum())

for name, n in [('control',161), ('physical_t2',802)]:
    d = root/name
    assert (d/'exit.code').read_text().strip() == '0'
    env = json.loads((d/'environment.json').read_text())
    assert env['EUNHA_KEPLER_SOFT'] == 'spline'
    for k,v in [('EUNHA_KEPLER_RC',.5),('EUNHA_KEPLER_CX',3),
                ('EUNHA_KEPLER_CY',3),('EUNHA_KEPLER_EVOLVE_R',2.5)]:
        assert float(env[k]) == v
    lines = (d/'log').read_text().splitlines()
    def records(tag):
        return [dict(x.split('=',1) for x in s.split()[1:])
                for s in lines if s.startswith('['+tag+']')]
    accepted = records('ACCEPTED_STATE')
    assert len(accepted)==1 and accepted[0]['finite']=='1' and accepted[0]['retries']=='0'
    dt = float(accepted[0]['dt'])
    stages = records('ORBIT_WORK_MEMORY_STAGE')
    legacy = records('ORBIT_WORK_STAGE')
    for rows in (stages,legacy):
        assert [s['stage'] for s in rows] == ['1','2','3','4']
        assert np.allclose([float(s['wdt']) for s in rows],
                           np.array([1,2,2,1])*dt/6, rtol=1e-14, atol=0)
    for s,l in zip(stages,legacy):
        assert s['wdt']==l['wdt'] and s['shear']==l['shear']
        vals = [float(s[k]) for k in ('hydro','shear','auxiliary','total')]
        assert np.isfinite(vals).all()
        assert np.isclose(sum(vals[:3]), vals[3], rtol=1e-13, atol=1e-30)
        assert np.isclose(float(l['hydro_minus_TCF'])+float(l['TCF']), vals[0], rtol=1e-13, atol=1e-30)
    integral = {k:sum(float(s['wdt'])*float(s[k]) for s in stages)
                for k in ('hydro','shear','auxiliary','total')}
    integral.update({k:sum(float(s['wdt'])*float(s[k]) for s in legacy)
                     for k in ('hydro_minus_TCF','TCF')})
    start = d/f'khout.{n:06d}.dat'; end = d/f'khout.{n+1:06d}.dat'
    a,t0 = read(start); c,t1 = read(end)
    assert np.isclose(t1-t0,dt,rtol=1e-11,atol=0)
    assert np.array_equal(a['id'],c['id']) and np.array_equal(a['mass'],c['mass'])
    ids0,e0 = departure(a); ids1,e1 = departure(c)
    assert np.array_equal(ids0,ids1)
    result = {'start_sha256':sha(start), 'end_sha256':sha(end),
              'start_time':t0, 'end_time':t1, 'RK_integrals':integral,
              'departure_change':e1-e0,
              'HG_plus_RK_remainder':e1-e0-integral['total']}
    if name == 'control':
        ref = root.parent/'Codex_hydro_partition_control_20261018'/end.name
        r,tr = read(ref)
        assert tr==t1 and all(np.array_equal(r[k],c[k]) for k in r.dtype.names)
        assert Path(str(ref)+'.rmem').read_bytes()==Path(str(end)+'.rmem').read_bytes()
        result['state_control_passed'] = True
    results[name] = result
results['scope'] = 'Completed one-step RK work analysis. Remainder combines finite HG split and RK update/quadrature error; not isolated HG work or long-term stability.'
print(json.dumps(results, indent=2))

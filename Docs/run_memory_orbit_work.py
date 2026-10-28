"""Same-stage memory/hydro departure-work diagnostic; two one-step branches.

The first reuses a completed reference to check diagnostic passivity. The
second reads a nonzero-memory t2 checkpoint. Neither changes the live run.
"""
import json, os, shutil, subprocess, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read, sha

b = Path('/gpfs/kjhan/LagForce')
root = b/'Codex_memory_orbit_work_control_20261018'
binary = b/'Codex_memory_orbit_work_build_20261018/eunha2.exe'
assert (binary.parent/'build.exit').read_text().strip() == '0'
spec = [('control', b/'Codex_hydro_partition_control_20261018', 161),
        ('physical_t2', b/'Codex_reversible_memory_physical_to_t10_20261018/memory_from_t0', 802)]
estimate = sum(2*((s/f'khout.{n:06d}.dat').stat().st_size +
                  (s/f'khout.{n:06d}.dat.rmem').stat().st_size) for _, s, n in spec)
assert shutil.disk_usage(b).free > 2*estimate
root.mkdir(exist_ok=False)
(root/'manifest.json').write_text(json.dumps({'binary': str(binary), 'sha256': sha(binary),
    'expected_source_and_final_bytes': estimate,
    'scope': 'Two one-step diagnostic branches, no long-time validation.'}, indent=2))
results = {}
for name, src, n in spec:
    d = root/name
    d.mkdir()
    snap = src/f'khout.{n:06d}.dat'
    _, t = read(snap)
    for p in (snap, Path(str(snap)+'.rmem'), src/'params.dat'):
        shutil.copyfile(p, d/p.name)
    case = json.loads((src/'environment.json').read_text())
    case.update(HYDRO_TSTOP=str(t+.00125), EUNHA_DUMP_DT='1000000',
                EUNHA_DUMP_FINAL='1', GFS_ORBIT_WORK_AUDIT='1')
    (d/'environment.json').write_text(json.dumps(case, indent=2))
    env = {k:v for k,v in os.environ.items() if not k.startswith(
        ('GFS_', 'EUNHA_', 'SEDOV_', 'LAGMFM_', 'HYDRO_', 'OMP_'))}
    env.update(case)
    with (d/'log').open('x') as f:
        p = subprocess.Popen(['timeout','900','mpirun','-np','1',str(binary),'params.dat',str(n)],
                             cwd=d, env=env, stdout=f, stderr=subprocess.STDOUT)
        (d/'pid').write_text(str(p.pid))
        print(name, p.pid, flush=True)
        code = p.wait()
    (d/'exit.code').write_text(str(code))
    assert code == 0
    lines = (d/'log').read_text().splitlines()
    accepted = [s for s in lines if s.startswith('[ACCEPTED_STATE]')]
    assert len(accepted)==1 and ' finite=1 ' in accepted[0] and ' retries=0 ' in accepted[0]
    stages = [dict(x.split('=',1) for x in s.split()[1:]) for s in lines
              if s.startswith('[ORBIT_WORK_MEMORY_STAGE]')]
    assert [s['stage'] for s in stages] == ['1','2','3','4']
    assert all(np.isfinite(float(s[k])) for s in stages for k in ('wdt','hydro','shear','auxiliary','total'))
    result = {'source':str(snap), 'source_sha256':sha(snap), 'start_time':t,
              'stages':stages, 'accepted':accepted[0],
              'RK_integrals':{k:sum(float(s['wdt'])*float(s[k]) for s in stages)
                              for k in ('hydro','shear','auxiliary','total')}}
    if name == 'control':
        a, ta = read(src/f'khout.{n+1:06d}.dat')
        c, tc = read(d/f'khout.{n+1:06d}.dat')
        assert ta == tc
        result['fields_bitwise_equal'] = {k:bool(np.array_equal(a[k],c[k])) for k in a.dtype.names}
        result['sidecar_bytes_equal'] = (src/f'khout.{n+1:06d}.dat.rmem').read_bytes() == (d/f'khout.{n+1:06d}.dat.rmem').read_bytes()
        assert all(result['fields_bitwise_equal'].values()) and result['sidecar_bytes_equal']
    results[name] = result
    (root/'analysis.json').write_text(json.dumps(results, indent=2))
    print(json.dumps(result, indent=2), flush=True)

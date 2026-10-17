"""Refine a completed 17-node Kepler audit, preserving all sampled nodes."""
import argparse,hashlib,json,os
from pathlib import Path
import numpy as np
from check_kepler_dense_step import atomic,simpson
from physical_field_accuracy_audit import environment
from check_kepler_moment_preflight import install
from adaptive_coupled_mass import restore
from gravity_coupled_2d import hydro_dense_state
from smooth_mass_path import evaluate_smooth_path
from kepler_disk_reference import frozen


def interleave_seed(rows,intervals):
    if len(rows)!=17 or intervals<32 or intervals%16:
        raise ValueError('require 17 seed nodes and an interval count divisible by 16')
    return {str(i*(intervals//16)):row for i,row in enumerate(rows)}


def run(seed,intervals=32):
    seed=Path(seed);data={n:(seed/n).read_bytes() for n in
        ('identity.json','trial.json','nodes.json','source.json','configuration.json')}
    old=json.loads(data['identity.json']);trial=json.loads(data['trial.json'])
    assert old['source_sha256']==hashlib.sha256(data['source.json']).hexdigest()
    assert old['configuration_sha256']==hashlib.sha256(data['configuration.json']).hexdigest()
    # Refine with the same physical modules; only this orchestration is new.
    def manifest(path):
        return {line.split(None,1)[1].strip():line.split()[0]
                for line in path.read_text().splitlines()}
    previous=manifest(seed.parent/'SHA256SUMS')
    current=manifest(Path(__file__).with_name('SHA256SUMS'))
    for name,digest in previous.items():
        if current.get(name)!=digest:raise ValueError('seed source changed: '+name)
    identity=dict(seed_files={n:hashlib.sha256(v).hexdigest() for n,v in data.items()},
        intervals=intervals,source_manifest_sha256=hashlib.sha256(Path(__file__).with_name('SHA256SUMS').read_bytes()).hexdigest())
    if Path('identity.json').exists():
        if json.loads(Path('identity.json').read_text())!=identity:raise ValueError('refinement identity mismatch')
    else:atomic('identity.json',identity)
    rows=json.loads(Path('nodes.json').read_text()) if Path('nodes.json').exists() else interleave_seed(json.loads(data['nodes.json']),intervals)
    os.environ.update(environment(json.loads(data['configuration.json'])))
    os.environ.update(LAGFORCE_VECTOR_KINETIC_ROOT='1',LAGFORCE_KINETIC_SMOOTHING='1',
        LAGFORCE_KINETIC_SMOOTHING_LENGTH='1',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star')
    install();record=restore(trial['record']);mass=np.asarray(json.loads(data['source.json'])['mass'])
    mask=np.asarray(record['frozen_first'],dtype=bool)
    for node in range(intervals+1):
        if str(node) in rows:continue
        s=node/intervals;state,F,v=hydro_dense_state(record,mass,s)
        if not np.array_equal(frozen(record['t']+s*record['dt'],state[:-4].reshape(-1,20)[:,:2]),mask):
            raise RuntimeError('frozen boundary crossing')
        rates=v[:-4].reshape(-1,20)
        flux,d=evaluate_smooth_path(record['t']+s*record['dt'],state,F,mass,
            (rates[:,:2],rates[:,2]),order=4,frozen_mask=mask)
        rows[str(node)]=dict(flux=flux,diag=d);atomic('nodes.json',rows)
        print('refined node',node,'/',intervals,flush=True)
    f=np.asarray([rows[str(i)]['flux'] for i in range(intervals+1)]);dt=record['dt'];active=~mask
    fine=simpson(f,dt);coarse=simpson(f[::2],dt)
    report=dict(identity=identity,dt=dt,reused_seed_nodes=17,
        absolute_relative_defect=float(np.max(simpson(abs(f),dt)[active]/mass[active])),
        signed_relative_defect=float(np.max(abs(fine[active])/mass[active])),
        nested_difference=float(np.max(abs(fine-coarse)[active]/mass[active])),
        maximum_relative_rate=float(np.max(abs(f[:,active])/mass[active])),
        full_validation=False)
    atomic('report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('seed');p.add_argument('--intervals',type=int,default=32)
    a=p.parse_args();run(a.seed,a.intervals)

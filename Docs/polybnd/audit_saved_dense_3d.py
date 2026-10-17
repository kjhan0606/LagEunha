"""Checkpointed quadrature of one saved actual 3D trajectory.

Prescribed geometry derivatives are never replaced by a fresh mass root.
Reports numerical estimates, not a continuum enclosure or whole-run pass.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from dense_record_3d import load
from moment_evolution_3d import dense_state
from kinetic_moment_rhs_3d import rhs
from power_cells_3d import Geometry3D


def atomic(path, data):
    path = Path(path)
    temp = Path(str(path) + '.tmp')
    temp.write_text(json.dumps(data, indent=2))
    temp.replace(path)


def integrate(values, dt):
    values = np.asarray(values)
    if len(values) < 3 or len(values) % 2 != 1:
        raise ValueError('Simpson quadrature requires an odd node count >= 3')
    return dt / (3 * (len(values)-1)) * (
        values[0]+values[-1]+4*values[1:-1:2].sum(axis=0)
        +2*values[2:-1:2].sum(axis=0))


def run(library, path, intervals, face_order=3):
    if intervals < 4 or intervals % 4:
        raise ValueError('intervals must be a positive multiple of four')
    record, endpoint = load(path)
    mass = record['raw'][:, 0]
    identity = dict(record_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                    intervals=intervals, face_order=face_order,
                    library_sha256=hashlib.sha256(Path(library).read_bytes()).hexdigest(),
                    source_sha256=hashlib.sha256(Path(__file__).with_name('SHA256SUMS').read_bytes()).hexdigest())
    saved = dict(identity=identity, rows=[])
    if Path('progress.json').exists():
        saved = json.loads(Path('progress.json').read_text())
        if saved['identity'] != identity:
            raise ValueError('audit identity mismatch')
    end = dense_state(record, 1.)
    for actual, key in zip(end[:4], ('sites', 'weights', 'raw', 'factors')):
        np.testing.assert_allclose(actual, endpoint[key], rtol=0, atol=3e-14)
    geo = Geometry3D(library)
    for node in range(len(saved['rows']), intervals+1):
        x,w,r,F,v,z = dense_state(record, node/intervals)
        np.testing.assert_array_equal(r[:, 0], mass)
        _,_,_,d = rhs(geo,x,w,r,F,np.array([[0.,0.,0.],[1.,1.,1.]]),
                       lambda p:np.zeros(3),face_order=face_order,pressure_reconstruction=record.get('pressure_reconstruction','quadratic'),face_integration=record.get('face_integration','point'),
                       audit_geometry_rates=(v,z))
        if d['mass_root_solved']:
            raise RuntimeError('diagnostic unexpectedly solved a mass root')
        saved['rows'].append(dict(rate=d['mass_rates'].tolist(),
                                  generator=d['generator_velocity_defect']))
        atomic('progress.json', saved)
        print('dense audit node',node,'/',intervals,flush=True)
    rates=np.asarray([r['rate'] for r in saved['rows']]);dt=record['dt']
    fine=integrate(rates,dt);coarse=integrate(rates[::2],dt)
    report=dict(identity=identity,t=record['t'],dt=dt,
                integrated_relative_defect=float(np.max(abs(fine)/mass)),
                integrated_absolute_relative_defect=float(np.max(integrate(abs(rates),dt)/mass)),
                nested_quadrature_difference=float(np.max(abs(fine-coarse)/mass)),
                maximum_relative_mass_rate=float(np.max(abs(rates)/mass)),
                generator_path_defect=max(r['generator'] for r in saved['rows']),
                whole_run_certified=False)
    atomic('report.json',report)
    return report


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('library');p.add_argument('record')
    p.add_argument('--intervals',type=int,default=16)
    p.add_argument('--face-order',type=int,default=3);a=p.parse_args()
    print(json.dumps(run(a.library,a.record,a.intervals,a.face_order)),flush=True)

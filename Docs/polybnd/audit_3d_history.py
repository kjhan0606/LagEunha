"""Audit every accepted saved 3D path up to a captured checkpoint.

Independent split mass quadrature; never solves for alternative weight rates.
Each sampled node is cached with its record/library/source identity. A partial
run is reported as a prefix, never as completion of the intended interval.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from dense_record_3d import load
from moment_evolution_3d import dense_state
from power_cells_3d import Geometry3D
from split_face_mass_3d import evaluate


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def atomic(path,data):
    path=Path(path);temp=Path(str(path)+'.tmp')
    temp.write_text(json.dumps(data,default=lambda x:x.tolist()));temp.replace(path)


def simpson(values,dt):
    return dt/(3*(len(values)-1))*(values[0]+values[-1]+4*values[1:-1:2].sum(0)+2*values[2:-1:2].sum(0))


def validate_chain(record,endpoint,previous,t,mass):
    if record['t']!=t:raise ValueError('dense path time gap')
    for k in ('sites','weights','raw','factors'):
        np.testing.assert_array_equal(record[k],previous[k])
    np.testing.assert_array_equal(record['raw'][:,0],mass)
    np.testing.assert_array_equal(endpoint['raw'][:,0],mass)


def summarize(steps,mass,duration):
    signed=np.zeros(len(mass));absolute=signed.copy();quad=signed.copy();peak=0.
    for row in steps:
        signed+=row['signed'];absolute+=row['absolute'];quad+=row['quad_difference']
        peak=max(peak,float(np.max(abs(signed)/mass)))
    return dict(max_absolute_relative_mass=float(np.max(absolute/mass)),
        max_signed_relative_mass=float(np.max(abs(signed)/mass)),
        max_prefix_signed_relative_mass=peak,
        summed_nested_quadrature_difference=float(np.max(quad/mass)),
        duration_scaled_peak_rate=duration*max((r['relative_rate'] for r in steps),default=0.),
        maximum_generator_path_defect=max((r['generator'] for r in steps),default=0.))


def run(run_dir,library,order=4,intervals=16):
    if intervals<4 or intervals%4:raise ValueError('interval count must be divisible by four')
    root=Path(run_dir);cfg=json.loads((root/'configuration.json').read_text())
    initial=json.loads((root/'initial_state.json').read_text())
    checkpoint_bytes=(root/'checkpoint.json').read_bytes();cp=json.loads(checkpoint_bytes)
    identity=dict(initial_sha256=digest(root/'initial_state.json'),configuration_sha256=digest(root/'configuration.json'),
        library_sha256=digest(library),source_manifest_sha256=digest(Path(__file__).with_name('SHA256SUMS')),order=order,intervals=intervals)
    if Path('identity.json').exists():
        if json.loads(Path('identity.json').read_text())!=identity:raise ValueError('history audit identity changed')
    else:atomic('identity.json',identity)
    mass=np.asarray(initial['mass']);previous={k:np.asarray(initial[k]) for k in ('sites','weights','raw','factors')}
    t=initial['t'];geo=Geometry3D(library);steps=[]
    for number in range(1,cp['steps']+1):
        p=root/f'step_{number:05d}.npz';record,endpoint=load(p)
        validate_chain(record,endpoint,previous,t,mass)
        cache=Path(f'step_{number:05d}_nodes.json');key=dict(identity=identity,record_sha256=digest(p))
        saved=json.loads(cache.read_text()) if cache.exists() else dict(key=key,rows=[])
        if saved['key']!=key:raise ValueError('node cache identity changed')
        for node in range(len(saved['rows']),intervals+1):
            x,w,a,L,v,z=dense_state(record,node/intervals)
            np.testing.assert_array_equal(a[:,0],mass)
            flux,d=evaluate(geo,x,w,a,L,(v,z),order=order,return_diagnostics=True)
            saved['rows'].append(dict(flux=flux,generator=d['generator_defect']));atomic(cache,saved)
            print('history audit step',number,'node',node,'/',intervals,flush=True)
        rates=np.asarray([r['flux'] for r in saved['rows']]);dt=record['dt'];signed=simpson(rates,dt)
        row=dict(signed=signed,absolute=simpson(abs(rates),dt),quad_difference=abs(signed-simpson(rates[::2],dt)),
            relative_rate=float(np.max(abs(rates)/mass)),generator=max(r['generator'] for r in saved['rows']))
        steps.append(row);previous=endpoint;t=record['t']+dt
        atomic('prefix.json',dict(t=t,steps=number,**summarize(steps,mass,cfg['end']-initial['t'])))
    if t!=cp['t']:raise ValueError('captured checkpoint time does not match path')
    for k in ('sites','weights','raw','factors'):np.testing.assert_array_equal(previous[k],cp[k])
    report=dict(t=t,target=cfg['end'],steps=cp['steps'],covers_target=bool(t>=cfg['end']-1e-15),
        captured_checkpoint_sha256=hashlib.sha256(checkpoint_bytes).hexdigest(),
        **summarize(steps,mass,cfg['end']-initial['t']),full_validation=False,
        note='Nested numerical quadrature estimates, not continuum enclosures; physical errors and other gates remain separate.')
    atomic('report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run_dir');p.add_argument('library')
    p.add_argument('--order',type=int,default=4);p.add_argument('--intervals',type=int,default=16)
    a=p.parse_args();run(a.run_dir,a.library,a.order,a.intervals)

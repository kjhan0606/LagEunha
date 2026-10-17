"""Audit one recorded split step from a preserved Kepler checkpoint.

The gravity kicks leave geometry/density moments fixed. This audits the actual
hydro subflow between them, with frozen cells treated as an external reservoir.
The trial endpoint is not promoted to an accepted production checkpoint.
"""
import argparse,hashlib,json,os
from pathlib import Path
import numpy as np
from physical_field_accuracy_audit import environment
from check_kepler_moment_preflight import install
from kepler_disk_reference import frozen,CENTER
from kepler_spline_gravity import acceleration
from gravity_coupled_2d import step,hydro_dense_state
from adaptive_coupled_mass import restore
from smooth_mass_path import evaluate_smooth_path


def atomic(name,data):
    path=Path(name);temp=Path(str(path)+'.tmp')
    temp.write_text(json.dumps(data,default=lambda x:x.tolist()))
    temp.replace(path)


def simpson(values,dt):
    return dt/(3*(len(values)-1))*(values[0]+values[-1]
        +4*values[1:-1:2].sum(axis=0)+2*values[2:-1:2].sum(axis=0))


def run(source,configuration,dt):
    data=Path(source).read_bytes();cfgdata=Path(configuration).read_bytes()
    cp=json.loads(data);cfg=json.loads(cfgdata)
    os.environ.update(environment(cfg))
    os.environ.update(LAGFORCE_VECTOR_KINETIC_ROOT='1',LAGFORCE_KINETIC_SMOOTHING='1',
        LAGFORCE_KINETIC_SMOOTHING_LENGTH='1',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star')
    install()
    identity=dict(source_sha256=hashlib.sha256(data).hexdigest(),
        configuration_sha256=hashlib.sha256(cfgdata).hexdigest(),dt=dt,
        source_manifest_sha256=hashlib.sha256(Path(__file__).with_name('SHA256SUMS').read_bytes()).hexdigest())
    if Path('identity.json').exists():
        if json.loads(Path('identity.json').read_text())!=identity:raise ValueError('identity mismatch')
    else:atomic('identity.json',identity)
    mass=np.asarray(cp['mass']);y=np.asarray(cp['state']);L=[np.asarray(a) for a in cp['factors']]
    if not Path('trial.json').exists():
        out,F,d=step(cp['t'],y,mass,L,dt,lambda p:acceleration(p,center=CENTER,rc=.5),
                     frozen=frozen,return_record=True)
        if d['stage_mass_flux']>=2e-13 or d['generator_velocity_defect']!=0:
            raise RuntimeError('trial stage constraint failure')
        atomic('trial.json',dict(record=d.pop('hydro_record'),endpoint=dict(t=cp['t']+dt,
            state=out,mass=mass,factors=F),diagnostics=d))
    trial=json.loads(Path('trial.json').read_text());record=restore(trial['record'])
    rows=json.loads(Path('nodes.json').read_text()) if Path('nodes.json').exists() else []
    for node in range(len(rows),17):
        s=node/16;state,F,derivative=hydro_dense_state(record,mass,s)
        mask=np.asarray(record['frozen_first'],dtype=bool)
        if not np.array_equal(frozen(cp['t']+s*dt,state[:-4].reshape(-1,20)[:,:2]),mask):
            raise RuntimeError('frozen boundary crossing requires piecewise audit')
        rates=derivative[:-4].reshape(-1,20)
        flux,d=evaluate_smooth_path(cp['t']+s*dt,state,F,mass,
            (rates[:,:2],rates[:,2]),order=4,frozen_mask=mask)
        rows.append(dict(flux=flux,diag=d));atomic('nodes.json',rows)
        print('Kepler hydro path',node,'/16',flush=True)
    flux=np.asarray([a['flux'] for a in rows]);active=~np.asarray(record['frozen_first'],dtype=bool)
    fine=simpson(flux,dt);coarse=simpson(flux[::2],dt)
    report=dict(identity=identity,t=cp['t'],dt=dt,stage=trial['diagnostics'],
        integrated_active_absolute_relative_defect=float(np.max(simpson(abs(flux),dt)[active]/mass[active])),
        integrated_active_signed_relative_defect=float(np.max(abs(fine[active])/mass[active])),
        active_nested_quadrature_difference=float(np.max(abs(fine-coarse)[active]/mass[active])),
        maximum_active_relative_rate=float(np.max(abs(flux[:,active])/mass[active])),
        frozen_reservoir_signed_mass_exchange=float(fine[~active].sum()),
        maximum_generator_path_defect=max(a['diag']['generator_defect'] for a in rows),
        full_validation=False,endpoint_accepted=False)
    atomic('report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('configuration')
    p.add_argument('--dt',type=float,required=True);a=p.parse_args();run(a.source,a.configuration,a.dt)

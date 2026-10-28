"""Checkpointed enriched Kepler continuation with a cached audit of every step.

Each fluid trial is saved before path quadrature. Resume reuses that trial and
all completed audit nodes. Failed gates leave the accepted checkpoint intact.
This auxiliary continuation does not certify the canonical 100-orbit problem.
"""
import argparse,hashlib,json,os
from pathlib import Path
import numpy as np
from audit_kepler_energy_step import configure,phi,atomic
from potential_moment_coupled_2d import step
from kepler_disk_reference import polar_rule,frozen,CENTER
from kepler_spline_gravity import acceleration
from audit_enriched_kepler_path import run as audit_path
from check_kepler_dense_step import simpson


def validation_duration(start,end,validation_end=None):
    horizon=end if validation_end is None else validation_end
    if not all(np.isfinite(x) for x in (start,end,horizon)) or not start<end or horizon<end:
        raise ValueError('finite validation horizon must cover execution interval')
    return horizon-start


def run(source,configuration,dt,end,gravity_scheme='finite',validation_end=None):
    source=Path(source).resolve();configuration=Path(configuration).resolve();folder=Path.cwd()
    data=source.read_bytes();cfg=configuration.read_bytes();initial=json.loads(data)
    duration=validation_duration(initial['t'],end,validation_end)
    identity=dict(source_sha256=hashlib.sha256(data).hexdigest(),configuration_sha256=hashlib.sha256(cfg).hexdigest(),
        source_manifest_sha256=hashlib.sha256(Path(__file__).with_name('SHA256SUMS').read_bytes()).hexdigest(),dt=dt,end=end,gravity_scheme=gravity_scheme)
    if validation_end is not None:identity['validation_end']=validation_end
    if Path('identity.json').exists():
        if json.loads(Path('identity.json').read_text())!=identity:raise ValueError('continuation identity changed')
    else:atomic('identity.json',identity)
    if 'potential_moments' not in initial:raise ValueError('continuation requires existing potential moments')
    if not np.isfinite(dt) or not 0<dt:raise ValueError('invalid timestep')
    configure(json.loads(cfg));cp=json.loads(Path('checkpoint.json').read_text()) if Path('checkpoint.json').exists() else initial
    mass=np.asarray(initial['mass']);a0=np.asarray(initial['state'])[:-4].reshape(-1,20)
    mask=frozen(initial['t'],a0[:,:2]);active=~mask
    thermal=float((a0[active,19]-.5*(a0[active,15]+a0[active,18])).sum())
    count=cp.get('accepted_steps',0);absolute=np.asarray(cp.get('path_absolute_cell',np.zeros(len(mass))))
    signed=np.asarray(cp.get('path_signed_cell',np.zeros(len(mass))));quad=np.asarray(cp.get('path_quad_cell',np.zeros(len(mass))))
    exchange=cp.get('potential_exchange',0.);peak_rate=cp.get('peak_relative_mass_rate',0.)
    while cp['t']<end-1e-15:
        np.testing.assert_array_equal(cp['mass'],mass)
        t=cp['t'];h=min(dt,end-t);y=np.asarray(cp['state']);L=[np.asarray(x) for x in cp['factors']];U=np.asarray(cp['potential_moments'])
        sub=folder/f'step_{count+1:05d}';sub.mkdir(exist_ok=True);trial_path=sub/'trial.json'
        start_hash=hashlib.sha256(json.dumps(cp,sort_keys=True).encode()).hexdigest()
        if not trial_path.exists():
            from coupled_stage_cache import StageCache
            cache=StageCache(sub/'stages',dict(run=identity,start_hash=start_hash,t=t,dt=h))
            out,F,newU,d=step(t,y,mass,L,U,h,phi,lambda x:acceleration(x,center=CENTER,rc=.5),
                lambda p:polar_rule(p,16,radii=(.25,.5)),frozen,return_record=True,gravity_scheme=gravity_scheme,stage_cache=cache)
            record=d.pop('hydro_record')
            atomic(trial_path,dict(t=t+h,state=out,mass=mass,factors=F,potential_moments=newU,
                diagnostics=d,hydro_record=record,start_hash=start_hash))
        trial=json.loads(trial_path.read_text())
        if trial.get('start_hash')!=start_hash:raise ValueError('pending trial source differs')
        audit_dir=sub/'path';audit_dir.mkdir(exist_ok=True)
        try:
            os.chdir(audit_dir);audit_path(trial_path,configuration)
        finally:os.chdir(folder)
        report=json.loads((audit_dir/'report.json').read_text());nodes=json.loads((audit_dir/'nodes.json').read_text())
        rates=np.asarray([x['flux'] for x in nodes]);inc=simpson(rates,h)
        new_abs=absolute+simpson(abs(rates),h);new_signed=signed+inc;new_quad=quad+abs(inc-simpson(rates[::2],h))
        new_peak=max(peak_rate,float(np.max(abs(rates[:,active])/mass[active])))
        new_exchange=exchange+report['integrated_potential_boundary']
        a=np.asarray(trial['state'])[:-4].reshape(-1,20);newU=np.asarray(trial['potential_moments'])
        np.testing.assert_array_equal(a[mask],a0[mask]);np.testing.assert_array_equal(newU[mask],np.asarray(initial['potential_moments'])[mask])
        for i in np.flatnonzero(mask):np.testing.assert_array_equal(trial['factors'][i],initial['factors'][i])
        gas_boundary=float(trial['state'][-2]-initial['state'][-2])
        energy=float((a[active,19]-a0[active,19]).sum()+(newU[active]-np.asarray(initial['potential_moments'])[active]).sum())-gas_boundary-new_exchange
        metrics=dict(absolute_mass=float(np.max(new_abs[active]/mass[active])),
            quadrature_mass=float(np.max(new_quad[active]/mass[active])),rate=duration*new_peak,
            generator=report['maximum_generator_defect'],energy_over_initial_thermal=abs(energy)/thermal)
        limits=dict(absolute_mass=1e-5,quadrature_mass=1e-7,rate=1e-4,generator=1e-5,energy_over_initial_thermal=1e-5)
        passed=all(metrics[k]<=limits[k] for k in limits) and trial['diagnostics']['stage_mass_flux']<2e-13 and trial['diagnostics']['generator_velocity_defect']==0
        atomic(sub/'acceptance.json',dict(metrics=metrics,limits=limits,passed=passed,validation_duration=duration,full_validation=False))
        if not passed:raise RuntimeError('enriched continuation gate failed; accepted checkpoint preserved')
        count+=1;absolute,signed,quad=new_abs,new_signed,new_quad;exchange=new_exchange;peak_rate=new_peak
        cp={k:trial[k] for k in ('t','state','mass','factors','potential_moments')}
        cp.update(accepted_steps=count,path_absolute_cell=absolute.tolist(),path_signed_cell=signed.tolist(),path_quad_cell=quad.tolist(),
            potential_exchange=exchange,peak_relative_mass_rate=peak_rate,last_metrics=metrics)
        atomic('checkpoint.json',cp);print('enriched accepted',count,cp['t'],metrics,flush=True)
    atomic('report.json',dict(completed=True,t=cp['t'],source_t=initial['t'],steps=count,full_validation=False,
        validation_end=initial['t']+duration,covers_validation_end=bool(cp['t']>=initial['t']+duration-1e-15),
        scope='audited auxiliary continuation; field accuracy and canonical longrun remain separate'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('configuration');p.add_argument('--dt',type=float,required=True);p.add_argument('--end',type=float,required=True);p.add_argument('--gravity-scheme',choices=('finite','midpoint'),default='finite')
    p.add_argument('--validation-end',type=float)
    a=p.parse_args();run(a.source,a.configuration,a.dt,a.end,a.gravity_scheme,a.validation_end)

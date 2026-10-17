"""Continue a saved coupled Kepler fluid state; no restart of its prefix.

Stage constraints are checked. Dense path and full potential-energy boundary
balance remain separate gates, explicitly not certified here.
"""
import argparse,json,os,hashlib,time,traceback
from pathlib import Path
import numpy as np
from physical_field_accuracy_audit import environment
from check_kepler_moment_preflight import install
from kepler_disk_reference import frozen,CENTER
from kepler_spline_gravity import acceleration
from gravity_coupled_2d import step


def run(source,configuration,dt,end):
    cfg=json.loads(Path(configuration).read_text());os.environ.update(environment(cfg));os.environ.update(LAGFORCE_VECTOR_KINETIC_ROOT='1',LAGFORCE_KINETIC_SMOOTHING='1',LAGFORCE_KINETIC_SMOOTHING_LENGTH='1',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star');install()
    identity=dict(source_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest(),configuration_sha256=hashlib.sha256(Path(configuration).read_bytes()).hexdigest(),dt=dt,end=end)
    if Path('configuration_receipt.json').exists():
        if json.loads(Path('configuration_receipt.json').read_text())!=identity:raise ValueError('continuation identity mismatch')
    else:Path('configuration_receipt.json').write_text(json.dumps(identity,indent=2))
    initial=json.loads(Path(source).read_text());cp=json.loads(Path('checkpoint.json').read_text()) if Path('checkpoint.json').exists() else initial
    mass=np.asarray(cp['mass']);np.testing.assert_array_equal(mass,initial['mass']);y=np.asarray(cp['state']);L=[np.asarray(a) for a in cp['factors']];t=cp['t'];count=cp.get('continuation_steps',0);peak=cp.get('stage_mass_peak',0.);start=time.monotonic()
    if not 0<dt or not initial['t']<=t<=end:raise ValueError('invalid continuation time')
    try:
        while t<end-1e-15:
            h=min(dt,end-t);mask=frozen(t,y[:-4].reshape(-1,20)[:,:2]);old=y.copy();oldL=[a.copy() for a in L]
            out,newL,d=step(t,y,mass,L,h,lambda p:acceleration(p,center=CENTER,rc=.5),frozen=frozen)
            np.testing.assert_array_equal(out[:-4].reshape(-1,20)[mask],old[:-4].reshape(-1,20)[mask])
            for i in np.flatnonzero(mask):np.testing.assert_array_equal(newL[i],oldL[i])
            if d['stage_mass_flux']>=2e-13 or d['generator_velocity_defect']!=0:raise RuntimeError('stage mass or steering gate failed')
            y,L=out,newL;t+=h;count+=1;peak=max(peak,d['stage_mass_flux'])
            cp=dict(t=t,state=y.tolist(),mass=mass.tolist(),factors=[a.tolist() for a in L],continuation_steps=count,stage_mass_peak=peak,last_diagnostics=d)
            Path('checkpoint.tmp').write_text(json.dumps(cp));Path('checkpoint.tmp').replace('checkpoint.json');print('Kepler continuation',count,t,flush=True)
    except Exception as exc:
        Path('failure.json').write_text(json.dumps(dict(t=t,error=repr(exc),traceback=traceback.format_exc()),indent=2));raise
    result=dict(completed=True,source_time=initial['t'],t=t,dt=dt,steps=count,stage_mass_peak=peak,dense_path_certified=False,full_validation=False,seconds=time.monotonic()-start)
    Path('report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('configuration');p.add_argument('--dt',type=float,required=True);p.add_argument('--end',type=float,required=True);a=p.parse_args();run(a.source,a.configuration,a.dt,a.end)

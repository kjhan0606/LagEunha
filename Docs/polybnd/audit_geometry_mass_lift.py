"""Improve preserved failed geometry paths without repeating approved evolution."""
import json,time
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic,restore
from saved_curved_contact import checked_json
from geometry_mass_lift import integrate,lifted_state,raw_path
from smooth_mass_path import sample_record
from smooth_path_quadrature import assess_values
from simple_wave_evolution import install_domain


def run():
    install_domain();cfg=json.loads(Path('geometry_lift_audit_config.json').read_text());rows=[]
    for item in cfg['cases']:
        started=time.monotonic();source=checked_json(item['initial'],item['initial_sha256']);mass=np.asarray(source['mass'])
        trial=checked_json(item['trial'],item['trial_sha256']);record=restore(trial['record'])
        record['geometry_lift']=integrate(record,mass)
        atomic('n'+str(item['n'])+'_lifted_record.json',dict(record=record))
        values=[sample_record(record,mass,j/16) for j in range(17)]
        old=[sample_record(trial['record'],mass,j/8) for j in range(9)]
        checked=assess_values(values,mass,record['dt'],.05,cfg['budgets'],True)
        previous=assess_values(old,mass,record['dt'],.05,cfg['budgets'],True)
        state,L,_=lifted_state(record,mass,1.);unchanged,_,_=raw_path(record,mass,1.);a=state[:-4].reshape(-1,20);oldend=np.asarray(trial['new'])[:-4].reshape(-1,20)
        result=dict(n=item['n'],dt=record['dt'],old_ratio=previous['ratio'],ratio=checked['ratio'],physical_ratio=checked['physical_ratio'],quadrature_ratio=checked['quadrature_ratio'],
            mass_ratio=checked['mass_ratio'],generator_defect=checked['generator'],maximum_absolute_mass_defect=float(np.max(checked['absolute']/mass)),
            root_peak=record['geometry_lift']['maximum_mass_flux'],lift_evaluations=record['geometry_lift']['evaluations'],
            endpoint_geometry_change=float(np.max(abs(a[:,:3]-oldend[:,:3]))),
            all_raw_path_moments_unchanged=bool(np.array_equal(a[:,3:],unchanged[:-4].reshape(-1,20)[:,3:])),
            endpoint_raw_roundoff_difference=float(np.max(abs(a[:,3:]-oldend[:,3:]))),seconds=time.monotonic()-started)
        rows.append(result);atomic('n'+str(item['n'])+'_lift_audit.json',result);print(result,flush=True)
    atomic('geometry_lift_audit_report.json',dict(completed=True,rows=rows,
        note='Only preserved rejected paths modified. Raw RK2 moments, positive factors and original masses unchanged. Independent saved-polynomial derivative audit, not a replacement root solve used as evidence.'))


if __name__=='__main__':run()

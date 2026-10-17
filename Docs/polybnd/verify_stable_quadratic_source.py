"""Revalidate the corrected numerical budget and preserved failed paths."""
import hashlib,json,multiprocessing as mp
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from adaptive_coupled_mass import atomic
from saved_curved_contact import checked_json
from simple_wave_evolution import install_domain
from smooth_mass_path import initialize_worker,sample_record
from smooth_path_quadrature import assess_values
import kinetic_factor_rhs as kinetic


def evaluate_task(task):
    record,mass,s=task
    return sample_record(record,np.asarray(mass),s)


def run():
    cfg=json.loads(Path('stable_quadratic_config.json').read_text());install_domain()
    source=checked_json(cfg['initial'],cfg['initial_sha256']);old=checked_json(cfg['rhs'],cfg['rhs_sha256'])
    mass=np.asarray(source['mass']);state=np.asarray(source['state']);L=[np.asarray(v) for v in source['factors']]
    dy,diag=kinetic.rhs(0.,state,mass,covariance_factors=L,face_order=4)
    difference=float(np.max(abs(dy-np.asarray(old['dy']))))
    atomic('stable_source_rhs.json',dict(t=0.,state_sha256=hashlib.sha256(state.tobytes()).hexdigest(),
        mass_sha256=hashlib.sha256(mass.tobytes()).hexdigest(),dy=dy,diag=diag))
    rows=[]
    with ProcessPoolExecutor(max_workers=8,mp_context=mp.get_context('spawn'),initializer=initialize_worker) as executor:
        for entry in cfg['failed_trials']:
            trial=checked_json(entry['path'],entry['sha256']);record=trial['record'];dt=record['dt']
            values=list(executor.map(evaluate_task,[(record,mass,j/8) for j in range(9)]))
            checked=assess_values(values,mass,dt,.002,cfg['budgets'],False)
            rows.append(dict(dt=dt,ratio=checked['ratio'],physical_ratio=checked['physical_ratio'],quadrature_ratio=checked['quadrature_ratio'],
                maximum_split_error=max(v['diag']['velocity_quadratic_covariance_split_error'] for v in values),
                maximum_orthogonality_error=max(v['diag']['velocity_quadratic_orthogonality_error'] for v in values)))
    passed=difference<1e-10 and diag['mass_flux']<2e-13 and all(r['maximum_split_error']<1e-9 for r in rows)
    report=dict(passed=bool(passed),source_rhs_difference=difference,source_mass_flux=diag['mass_flux'],rows=rows,
        note='Stable bounded SVD construction; same mathematical closure and budgets. One approved source RHS revalidated; old failed dynamical stages only read for independent path tests, no accepted time integration repeated.')
    atomic('stable_quadratic_validation.json',report);print(report,flush=True)
    raise SystemExit(0 if passed else 2)


if __name__=='__main__':run()

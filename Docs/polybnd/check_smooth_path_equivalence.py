"""Reproduce immutable serial samples before parallel continuation."""
import hashlib,json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from smooth_mass_path import initialize_worker,sample_record_file
from predictive_evolution_batch import record_digest
from saved_curved_contact import checked_json
from smooth_path_quadrature import assess_values


def run():
    cfg=json.loads(Path('parallel_equivalence_config.json').read_text());initialize_worker();rows=[]
    for entry in cfg['cases']:
        state=checked_json(entry['initial'],entry['initial_sha256'])
        trial=checked_json(entry['trial'],entry['trial_sha256']);samples=checked_json(entry['samples'],entry['samples_sha256'])
        mass=np.asarray(state['mass']);digest=record_digest(trial['record'])
        for j in (0,3,8):
            value=sample_record_file(entry['trial'],digest,mass,j/8);old=samples[str(j)]
            flux=float(np.max(abs(value['flux']-old['flux'])))
            diag=max(abs(value['diag'][k]-old['diag'][k]) for k in value['diag'])
            rows.append(dict(n=entry['n'],sample=j,flux_difference=flux,diagnostic_difference=diag))
        if 'report' in entry:
            report=checked_json(entry['report'],entry['report_sha256'])
            budgets=checked_json(entry['config'],entry['config_sha256'])
            checked=assess_values([samples[str(j)] for j in range(9)],mass,trial['record']['dt'],budgets['end'],budgets,False)
            same=checked['ratio']==report['ratio'] and all(np.array_equal(checked[k],report[k]) for k in ('absolute','signed','quad_delta'))
            rows.append(dict(n=entry['n'],sample='legacy_assessment',flux_difference=0. if same else 1.,diagnostic_difference=0.))
    result=dict(passed=all(v['flux_difference']==0 and v['diagnostic_difference']==0 for v in rows),rows=rows,
        note='Only six immutable audit samples recomputed for exact execution equivalence; no dynamical stage or accepted integration repeated.')
    atomic('smooth_path_equivalence.json',result);print(result,flush=True)
    raise SystemExit(0 if result['passed'] else 2)


if __name__=='__main__':run()

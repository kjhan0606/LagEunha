"""Resume preserved approved/pending work without repeating its source RHS."""
import json,hashlib,os
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from predictive_evolution_batch import advance_interval
from continue_common_contact import score


def resume_task(case,source,checkpoint):
    cfg=checkpoint['config']
    if cfg['case']!=case or cfg['interval']!='short' or cfg['start']!=source['t']:
        raise ValueError('saved branch identity mismatch')
    digest=hashlib.sha256(np.asarray(source['state']).tobytes()).hexdigest()
    if cfg['source_state_sha256']!=digest:raise ValueError('saved branch source state mismatch')
    if not np.array_equal(source['mass'],checkpoint['mass']):raise ValueError('saved mass array mismatch')
    if not cfg['start']<=checkpoint['t']<=cfg['end'] or cfg['duration']<=0 or cfg['initial_steps']<1:
        raise ValueError('invalid saved interval')
    for key,flag,enabled,disabled in [
        ('cap_inverse','LAGFORCE_ANALYTIC_CAP','analytic','legacy_bisection'),
        ('moment_fit_frame','LAGFORCE_CENTERED_MOMENT_FIT','cell_centered','global_raw'),
        ('power_cell_backend','LAGFORCE_ACTIVE_POWER_CELLS','active_all_constraints','exhaustive')]:
        actual=enabled if os.environ.get(flag,'0')=='1' else disabled
        if cfg.get(key,disabled)!=actual:raise ValueError('saved numerical option mismatch: '+key)
    return case,'short',cfg['duration'],cfg['initial_steps']


def run():
    config=json.loads(Path('resume_branch_config.json').read_text());case=config['case']
    if '/' in case or case in ('.','..'):raise ValueError('invalid case')
    source=json.loads((Path('inputs')/case/'checkpoint.json').read_text())
    path=Path('evolution')/case/'short';checkpoint=json.loads((path/'checkpoint.json').read_text())
    task=resume_task(case,source,checkpoint)
    result=advance_interval(task)
    report=dict(completed=True,passed=bool(result['passed']),result=result,
        source_time=source['t'],resumed_time=checkpoint['t'],
        note=config['note'],source_rhs_repeated=False)
    if result['passed']:
        final=json.loads((path/'checkpoint.json').read_text());mass=np.asarray(source['mass'])
        np.testing.assert_array_equal(final['mass'],mass)
        report['final_score']=score(final,source['config'],mass,'final')
        report['passed']=bool(report['passed'] and report['final_score']['reference_valid'])
    atomic('resumed_branch_report.json',report);print('resumed branch',case,report['passed'],flush=True)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()

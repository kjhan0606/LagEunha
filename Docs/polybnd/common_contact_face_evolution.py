"""Change face-bin count in actual evolution from the same preserved state."""
import hashlib
import json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from predictive_evolution_batch import advance_interval,initialize_case,density_options
from common_contact_time_refinement import digest,difference,uniform_grid
from continue_common_contact import score


def run():
    cfg=json.loads(Path('face_evolution_config.json').read_text())
    for path,expected in cfg['source_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=expected:
            raise RuntimeError('face experiment source changed: '+path)
    source=json.loads(Path(cfg['source_checkpoint']).read_text())
    initial=json.loads(Path(cfg['initial_config_checkpoint']).read_text())
    baseline=json.loads(Path(cfg['baseline_checkpoint']).read_text())
    baseline_result=json.loads(Path(cfg['baseline_result']).read_text())
    if not baseline_result['passed'] or baseline['config']['source_state_sha256']!=digest(source['state']):
        raise RuntimeError('baseline does not start at the same accepted source')
    original=initial['config'];mass,_=initialize_case(original)
    np.testing.assert_array_equal(source['mass'],mass)
    np.testing.assert_array_equal(baseline['mass'],mass)
    if baseline['config']['density_provider']!=density_options(original,True):
        raise RuntimeError('baseline reconstruction differs')
    changed=dict(original,density_provider=dict(original['density_provider'],segments=cfg['segments']))
    name='segments'+str(cfg['segments']);path=Path('inputs')/name/'checkpoint.json'
    if path.exists():raise RuntimeError('refuse duplicate face experiment input')
    path.parent.mkdir(parents=True)
    atomic(path,dict(source,config=changed,provenance=dict(source=cfg['source_checkpoint'],
        purpose='Only face segment count changes; actual continuation, not frozen-state scoring.')))
    result=advance_interval((name,'short',cfg['duration'],cfg['steps']))
    atomic('face_evolution_result.json',result)
    if not result['passed']:raise SystemExit(2)
    root=Path('evolution')/name/'short';final=json.loads((root/'checkpoint.json').read_text())
    np.testing.assert_array_equal(final['mass'],mass)
    if abs(final['t']-baseline['t'])>1e-15:raise RuntimeError('face endpoints differ in time')
    rows=[score(source,original,mass,'source64'),score(source,changed,mass,'source256'),
        score(baseline,original,mass,'evolved64'),score(final,changed,mass,'evolved256')]
    baseline_root=Path(cfg['baseline_checkpoint']).parent
    def grid(cp,folder):
        sizes=[json.loads((folder/('step_'+str(i)+'.json')).read_text())['dt'] for i in range(1,cp['steps']+1)]
        return uniform_grid(cp,sizes,cfg['duration'],cfg['steps'])
    same_grid=grid(baseline,baseline_root) and grid(final,root)
    report=dict(completed=True,passed=all(row['reference_valid'] for row in rows),
        source_time=source['t'],end_time=final['t'],segments=[64,cfg['segments']],
        actual_steps=[baseline['steps'],final['steps']],same_uniform_time_grid=same_grid,
        endpoint_differences=difference(baseline,final,mass),rows=rows,
        note='Same saved source, fixed actual mass, unsteered sites; only face segment count intentionally changes. If accepted time grids differ, this is not an isolated spatial discretization comparison. References used only in scoring.')
    atomic('face_evolution_report.json',report)
    print('face evolution comparison',report['passed'],'same time grid',same_grid,flush=True)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()

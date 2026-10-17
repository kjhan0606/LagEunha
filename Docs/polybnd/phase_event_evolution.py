"""Local numerical comparison through phase coalescence, checkpoint continuation."""
import json
from pathlib import Path
from adaptive_coupled_mass import atomic
from predictive_evolution_batch import advance_interval


def run():
    if not json.loads(Path('smoothing_report.json').read_text())['passed']:raise RuntimeError('smoothing validation failed')
    config=json.loads(Path('event_config.json').read_text());cp=json.loads(Path('inputs/n4_r100/checkpoint.json').read_text())
    if not cp['config']['smooth_phase_floor'] or cp['config']['path_control_horizon']!=1e-4:raise RuntimeError('event configuration mismatch')
    result=advance_interval(('n4_r100','short',config['target_time']-cp['t'],4))
    atomic('event_evolution_report.json',dict(passed=result['passed'],result=result,
        note='Separate local model comparison from a preserved pre-event state. Active unsmoothed AQ remains the reference. Same1e-4 path-control horizon and all original gates; source state/masses preserved.'))
    raise SystemExit(0 if result['passed'] else 2)


if __name__=='__main__':run()

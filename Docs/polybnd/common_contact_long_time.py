"""Longer fixed-model dt comparison from one preserved common-contact state."""
import hashlib
import json
from pathlib import Path
from adaptive_coupled_mass import atomic
from predictive_evolution_batch import advance_interval
import common_contact_time_refinement as comparison
from preserved_experiment_input import install_input


def run():
    config=json.loads(Path('long_time_config.json').read_text())
    for path,expected in config['source_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=expected:
            raise RuntimeError('preserved time source changed: '+path)
    source=json.loads(Path(config['source_checkpoint']).read_text())
    initial=json.loads(Path(config['initial_config_checkpoint']).read_text())
    cfg=dict(initial['config'],path_control_horizon=config['duration'])
    effective_initial=dict(initial,config=cfg)
    atomic('effective_initial_checkpoint.json',effective_initial)
    name='coarse'+str(config['coarse_steps'])
    path=Path('inputs')/name/'checkpoint.json'
    install_input(path,dict(source,config=cfg,provenance=dict(source=config['source_checkpoint'],
        purpose='New continuation from saved state; same physical model, longer horizon explicitly recorded.')))
    result=advance_interval((name,'short',config['duration'],config['coarse_steps']))
    atomic('coarse_evolution_report.json',result)
    if not result['passed']:raise SystemExit(2)
    root=Path('evolution')/name/'short'
    paths=[config['source_checkpoint'],'effective_initial_checkpoint.json',
        str(root/'checkpoint.json'),str(root/'result.json')]
    time_cfg=dict(source_checkpoint=paths[0],initial_config_checkpoint=paths[1],
        coarse_checkpoint=paths[2],coarse_result=paths[3],duration=config['duration'],
        fine_steps=config['fine_steps'],source_sha256={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths},
        purpose='Longer dt-only comparison; same preserved source/model/grid/face bins. Coarse result reused after its actual completion.')
    atomic('time_config.json',time_cfg)
    comparison.run()


if __name__=='__main__':run()

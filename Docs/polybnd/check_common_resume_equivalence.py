"""Verify saved path nodes still use the same equations before CPU migration."""
import json,time
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from predictive_evolution_batch import initialize_case,density_options,path_sample
from predictive_boundary_density import BoundaryDensity


def run():
    cfg=json.loads(Path('resume_probe_config.json').read_text());rows=[]
    for entry in cfg['entries']:
        case=entry['case'];cp=json.loads((Path('inputs')/case/'checkpoint.json').read_text())
        mass,_=initialize_case(cp['config']);data=json.loads(Path(entry['record_file']).read_text());record=data['record']
        nodes=data['samples'] if 'samples' in data else json.loads(Path(entry['nodes_file']).read_text())
        index=entry['index'];started=time.monotonic()
        actual=path_sample(record,mass,index,BoundaryDensity(**density_options(cp['config'],True)))
        error=float(np.max(abs(np.asarray(actual['mass_rate'])-nodes[str(index)]['mass_rate'])/mass))
        row=dict(case=case,index=index,max_relative_rate_difference=error,passed=error<1e-11,seconds=time.monotonic()-started)
        rows.append(row);print('resume node',case,index,error,flush=True)
    report=dict(passed=all(r['passed'] for r in rows),results=rows,note='Same equations, state and prescribed rates; only independent path-worker count changes. Saved nodes/accepted steps retained.')
    atomic('resume_equivalence_report.json',report);raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()

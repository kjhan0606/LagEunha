"""Extend completed coupled time convergence; do not rerun coarse cases."""
import json
from pathlib import Path
from coupled_factor_path_audit import run

if __name__=='__main__':
    prior=json.loads(Path('previous_summary.json').read_text())
    rows=prior['results']
    if [r['steps'] for r in rows]!=[16,32,64] or not all(r['status']=='completed' and r['state_gate'] for r in rows):
        raise RuntimeError('coarse coupled predecessor failed')
    results=[]
    for steps in (256,512,1024):
        r=run(steps);results.append(r)
        Path('refinement_summary.json').write_text(json.dumps(dict(previous=rows,results=results),indent=2))
        if r.get('status')!='completed' or not r.get('state_gate',False):break
        if r.get('mass_path_gate',False):break

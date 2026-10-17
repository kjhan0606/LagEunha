"""Analytic-Jacobian regression, then N3 constrained-path checks."""
import json
from pathlib import Path
import numpy as np
from instantaneous_weight_path import run
PREVIOUS=Path('/gpfs/kjhan/LagForce/Codex_G_instant_weights_20261010_081115')
def emit(r):
    with open('results.jsonl','a') as f:f.write(json.dumps(r)+'\n')
    print(json.dumps(r),flush=True)
def strict_predecessor_gate(rows):
    strict=[r for r in rows if r.get('rtol')==1e-10]
    expected={(kind,h) for kind in ('shear','wave') for h in (.01,.05)}
    coverage={(r.get('kind'),r.get('duration')) for r in strict}
    return len(strict)==4 and coverage==expected and all(r.get('status')=='completed' and r.get('sampled_path_gate') is True for r in strict)

def main():
    rows=[json.loads(s) for s in (PREVIOUS/'results.jsonl').read_text().splitlines()]
    if not strict_predecessor_gate(rows):
        Path('gate.json').write_text(json.dumps(dict(status='blocked_by_predecessor',records=rows),indent=2));return
    Path('predecessor_selection.json').write_text(json.dumps(dict(rtol=1e-10,reason='Proceed only with the four previously passed strict-tolerance cases; coarse failure remains recorded.',excluded=[r for r in rows if r.get('rtol')!=1e-10]),indent=2))
    r=run('shear',.01,1e-10,n=2,analytic=True);emit(r)
    old=json.loads((PREVIOUS/'shear_h0.01_tol1e-10.path.json').read_text())
    new=json.loads(Path('shear_n2_h0.01_tol1e-10_analytic1.path.json').read_text())
    errors={k:float(np.max(abs(np.array(old[k])-np.array(new[k])))) for k in ('final_flux','boundary')}
    errors['endpoint_weights']=float(np.max(abs(np.array(old['accepted_weights'][-1])-new['accepted_weights'][-1])))
    passed=r['sampled_path_gate'] and max(errors.values())<1e-10
    Path('jacobian_regression.json').write_text(json.dumps(dict(pass_gate=passed,differences=errors,old_calls=old['result']['calls'],new_calls=r['calls']),indent=2))
    if not passed:return
    for kind in ('shear','wave'):
        for h in (.01,.05):
            try:r=run(kind,h,1e-10,n=3,analytic=True)
            except Exception as e:r=dict(kind=kind,n=3,duration=h,status='failed',error=repr(e))
            emit(r)
if __name__=='__main__':main()

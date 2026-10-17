"""Instrument a preserved failed trial without changing its old source code."""
import importlib.util,json,sys,linecache,os
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from smooth_mass_path import sample_record,initialize_worker


def run():
    cfg=json.loads(Path('quadratic_budget_diagnosis_config.json').read_text());root=Path(cfg['source'])
    path=root/'quadratic_velocity_reconstruction.py'
    spec=importlib.util.spec_from_file_location('quadratic_velocity_reconstruction',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    sys.modules['quadratic_velocity_reconstruction']=module
    os.environ['LAGFORCE_QUADRATIC_VELOCITY']='1';initialize_worker()
    original=json.loads((root/'initial_state.json').read_text());trial=json.loads((root/'pending.json').read_text())
    mass=np.asarray(original['mass']);bad=[];failures=[]
    def trace(frame,event,arg):
        if frame.f_code.co_filename!=str(path) or frame.f_code.co_name!='prepare':return trace
        if event=='line' and 'orthogonal_error=max' in linecache.getline(str(path),frame.f_lineno):
            d=frame.f_locals
            error=float(np.linalg.norm(d['resolved']+d['L']@d['residual']@d['residual'].T@d['L'].T-d['R']))/d['norm']
            if error>1e-9:
                bad.append(dict(cell=d['i'],relative_split_error=error,covariance_norm=d['norm'],factor_condition=float(np.linalg.cond(d['L'])),candidate_singular_values=d['values']))
                if not Path('bad_budget_cell.npz').exists():
                    np.savez('bad_budget_cell.npz',L=d['L'],candidate=d['candidate'],gram_factor=d['gram_factor'],mass=d['m'])
        if event=='exception' and isinstance(arg[1],RuntimeError):
            d=frame.f_locals;failures.append(dict(split_error=d.get('split_error'),orthogonal_error=d.get('orthogonal_error'),message=str(arg[1])))
        return trace
    samples=[]
    for j in range(9):
        sys.settrace(trace)
        try:value=sample_record(trial['record'],mass,j/8);samples.append(dict(node=j,passed=True))
        except RuntimeError as exc:
            samples.append(dict(node=j,passed=False,error=str(exc)));break
        finally:sys.settrace(None)
    result=dict(source=str(root),samples=samples,failures=failures,bad_cells=bad,
        note='Preserved old code instrumented in an isolated process; old stages reused and no time integration repeated.')
    atomic('quadratic_budget_diagnosis.json',result);print(result,flush=True)


if __name__=='__main__':run()

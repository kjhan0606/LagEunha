"""Split accepted energy changes using existing stage logs, without rerunning."""
import argparse
import json
import math
import re
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('folder', type=Path)
a = p.parse_args()
rows = []
state = budget = None
def fields(line):
    return dict(re.findall(r'(\w+)=([^\s]+)', line))
with (a.folder / 'log').open() as f:
    for line in f:
        if not line.endswith('\n'):
            break
        if line.startswith('[ACCEPTED_STATE]'):
            assert state is None
            state = fields(line)
            assert state['finite'] == '1'
        elif line.startswith('[ACCEPTED_RK_BUDGET]'):
            budget = fields(line)
            assert state and budget['step'] == state['step']
            assert budget['stage_mask_min'] == budget['stage_mask_max'] == '15'
        elif line.startswith('[RMEM_ACCEPTED]'):
            z = fields(line)
            assert state and budget
            assert z['dt'] == state['dt'] and z['retries'] == state['retries']
            rk = float(budget['rk_dEh']) + float(budget['rk_dPhi'])
            hg = float(state['hg_dE'])
            delta = float(state['dEstep'])
            residual = float(budget['residual_E'])
            row = dict(step=int(state['step']), time=float(state['elapsed']),
                       hydro_change=delta, rk_stage_energy=rk, hg_change=hg,
                       rk_update_residual=residual,
                       identity_roundoff=delta-rk-hg-residual,
                       auxiliary_change=float(z['dEaux']),
                       auxiliary_boundary_work=float(z['boundary_work']))
            assert all(math.isfinite(v) for v in row.values())
            rows.append(row)
            state = budget = None
assert rows
keys = ('hydro_change', 'rk_stage_energy', 'hg_change', 'rk_update_residual',
        'identity_roundoff', 'auxiliary_change', 'auxiliary_boundary_work')
totals = {k: math.fsum(r[k] for r in rows) for k in keys}
totals['rk_plus_aux_minus_aux_boundary'] = (totals['rk_stage_energy']
    + totals['auxiliary_change'] - totals['auxiliary_boundary_work'])
print(json.dumps(dict(complete_steps=len(rows), latest=rows[-1], cumulative=totals,
    max_abs_step_rk_update_residual=max(abs(r['rk_update_residual']) for r in rows),
    largest_residual_step=max(rows, key=lambda r: abs(r['rk_update_residual'])),
    scope='Actual accepted RK4 increments include auxiliary acceleration. This '
    'separates update residual from RHS-integrated energy change, but does not '
    'separate hydro boundary work from spatial operator error. Not a conservation '
    'or stability pass. Complete state/RK/memory triples only.'), indent=2))

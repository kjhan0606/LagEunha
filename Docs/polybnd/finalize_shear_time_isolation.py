import json
from pathlib import Path
from adaptive_coupled_mass import atomic


def sensitivity(rows,baseline):
    values=[r['interface']['cell_average_density_L1_relative'] for r in rows]
    base=baseline['cell_average_density_L1_relative']
    if base<=0:raise ValueError('positive baseline error required')
    return dict(density_errors=values,adaptive_density_error=base,
        relative_changes_from_adaptive=[abs(v-base)/base for v in values],
        finest_change_relative_to_bias=abs(values[-1]-base)/base,
        note='A small finest change supports temporal insensitivity at this spatial resolution, not a proof of the remaining error source.')


def main():
    baseline=json.loads(Path('predecessor_completion.json').read_text())
    cases={}
    for n in (3,4):
        try:
            rows=[json.loads((Path(f'n{n}_{s}')/'result.json').read_text()) for s in (32,64,128)]
            if not all(r['interface']['reference_valid'] for r in rows):raise RuntimeError('reference invalid')
            cases[f'n{n}']=dict(status='completed',rows=rows,
                comparison=sensitivity(rows,baseline['cases'][f'n{n}_r100']))
        except Exception as exc:cases[f'n{n}']=dict(status='failed',error=repr(exc))
    result=dict(completed=all(c['status']=='completed' for c in cases.values()),cases=cases)
    atomic('completion_report.json',result);print(json.dumps(result),flush=True)
    raise SystemExit(0 if result['completed'] else 2)


if __name__=='__main__':main()

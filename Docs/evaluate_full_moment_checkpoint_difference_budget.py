"""Apply documented temporal or platform physical budgets to saved differences.

This never changes backend/seed roundoff bounds or authorizes an extension.
"""
import argparse,hashlib,json,math
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--comparison',type=Path,required=True)
    p.add_argument('--kind',choices=('platform','temporal'),required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();assert not args.output.exists()
    source=json.loads(args.comparison.read_text())
    assert source['status']=='completed_checkpoint_actual_cell_error_difference_diagnostics'
    assert source['existing_physical_monitors_both_passed'], 'Failed-monitor diagnostics cannot pass an acceptance budget'
    budget={'pressure':.0002,'density':.0002,'velocity':.0002,'phase_degrees':.02} if args.kind=='platform' else {'pressure':.002,'density':.002,'velocity':.002,'phase_degrees':.2}
    rows={}
    for name,row in source['regions'].items():
        phase=row['phase_error_difference_degrees']
        metrics={'pressure':row['relative_pressure_error_difference_RMS'],
            'density':row['relative_density_error_difference_RMS'],
            'velocity':row['velocity_error_difference_relative_RMS'],
            'phase_degrees':None if phase is None else phase['RMS']}
        passed={key:(value is None and key=='phase_degrees') or
            (value is not None and math.isfinite(value) and value<budget[key])
            for key,value in metrics.items()}
        rows[name]=dict(metrics=metrics,components_passed=passed,passed=all(passed.values()),
            cells_A=row['cells_A'],cells_B=row['cells_B'],intersection_cells=row['intersection_cells'],
            phase_excluded_cells=row['phase_excluded_cells'])
    assert rows and 'all' in rows
    report=dict(status=('passed_' if all(v['passed'] for v in rows.values()) else 'failed_')+args.kind+'_physical_difference_budget',
        time=source['time'],kind=args.kind,budget=budget,strict_less_than=True,regions=rows,
        comparison_path=str(args.comparison.resolve()),comparison_sha256=hashlib.sha256(args.comparison.read_bytes()).hexdigest(),
        automatic_extension_authorized=False,backend_or_seed_bounds_changed=False,
        physical_providers=0,accepted_time_steps=0,
        plan='Docs/astra_kepler_plan_review_2026-10-18.md density and platform budget amendment 2026-10-19',
        limits='Physical error difference only at this checkpoint. Missing times, individual density monitors, restart equivalence and contact growth still require separate evaluation.')
    args.output.mkdir();(args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],time=report['time'],failed_regions=[k for k,v in rows.items() if not v['passed']])),flush=True)


if __name__=='__main__':main()

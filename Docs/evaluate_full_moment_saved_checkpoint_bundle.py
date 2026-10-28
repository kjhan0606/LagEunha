"""Apply documented density/temporal budgets to saved checkpoints; no evolution.

This report never authorizes automatic extension. Missing checkpoints stay pending.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from compare_full_moment_checkpoint_errors import load, compare, REGIONS


def stats(error, mass, ids, mask):
    x=error[mask].astype(np.longdouble); w=mass[mask].astype(np.longdouble)
    j=int(np.argmax(abs(x)))
    return dict(cells=int(mask.sum()),RMS=float(np.sqrt(np.sum(w*x*x)/w.sum())),
                signed_mean=float(np.sum(w*x)/w.sum()),maximum=float(abs(x[j])),
                maximum_signed=float(x[j]),maximum_id=int(ids[mask][j]))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--initial-cavity',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    assert not args.output.exists()
    with np.load(args.initial_cavity) as f:
        initial_mask=f['initial_uniform_cavity'];initial_ids=f['ids']
    saved={};pending=[];pairs=[]
    for branch,steps in [('coarse',(2,5,10)),('fine',(4,10,20))]:
        rows=[]
        for step in steps:
            folder=args.root/branch/'checkpoints'/f'{branch}_{step:04d}'/'physics'
            # A completed physical gate is required; transient output is pending.
            if not (folder.parent/'physical_gate.json').exists():
                pending.append(str(folder));continue
            d,r,source=load(folder)
            assert np.array_equal(initial_ids,d['ids'])
            mass=d['actual_mass'];regions={}
            for name in REGIONS:
                mask=d[name]
                if not mask.any():continue
                density=stats(d['relative_density_reference_error'],mass,d['ids'],mask)
                regions[name]=dict(density=density,
                    density_passed=density['RMS']<.01 and density['maximum']<.1,
                    pressure=stats(d['relative_pressure_error'],mass,d['ids'],mask))
            masks={'initial_cavity':initial_mask,'current_cavity':d['cavity'],
                   'cavity_intersection':initial_mask&d['cavity'],
                   'left_cavity':initial_mask&~d['cavity'],
                   'entered_cavity':~initial_mask&d['cavity']}
            cavity={name:stats(d['relative_density_reference_error'],mass,d['ids'],mask)
                    if mask.any() else dict(cells=0) for name,mask in masks.items()}
            row=dict(time=r['time'],source=source,regions=regions,cavity=cavity,
                     density_passed=all(x['density_passed'] for x in regions.values()))
            rows.append(row);saved[(branch,r['time'])]=(d,row)
        # Signed/max-ID reports accompany growth ratios; ratio is no causal proof.
        for previous,current in zip(rows,rows[1:]):
            current['growth_since_previous']={'previous_time':previous['time'],'regions':{
                name:{field:current['regions'][name][field]['RMS']/previous['regions'][name][field]['RMS']
                      if previous['regions'][name][field]['RMS']>0 else None
                      for field in ('density','pressure')}
                for name in current['regions'] if name in previous['regions']}}
    for t in (.01,.025,.05):
        if ('coarse',t) not in saved or ('fine',t) not in saved:continue
        a,ra=saved[('coarse',t)];b,rb=saved[('fine',t)]
        regions=compare(a,b)
        for row in regions.values():
            phase=row['phase_error_difference_degrees']
            row['temporal_budget_passed']=(row['relative_pressure_error_difference_RMS']<.002
                and row['relative_density_error_difference_RMS']<.002
                and row['velocity_error_difference_relative_RMS']<.002
                and (phase is None or phase['RMS']<.2))
        pairs.append(dict(time=t,regions=regions,passed=all(x['temporal_budget_passed'] for x in regions.values())))
    report=dict(status='completed_saved_checkpoint_budget_assessment',
        plan='Docs/astra_kepler_plan_review_2026-10-18.md density/growth evaluation 2026-10-19',
        physical_density_budget=dict(RMS=.01,maximum=.1,strict_less_than=True),
        temporal_budget=dict(pressure_RMS=.002,density_RMS=.002,velocity_RMS=.002,phase_degrees_RMS=.2),
        checkpoints={branch:[row for (br,t),(d,row) in saved.items() if br==branch] for branch in ('coarse','fine')},
        comparisons=pairs,pending=pending,all_three_times_available=len(pairs)==3,
        available_density_passed=all(row['density_passed'] for d,row in saved.values()),
        available_temporal_passed=all(row['passed'] for row in pairs),
        automatic_extension_authorized=False,
        extension_reason='Completed t0.05 physical, temporal, platform and restart evidence plus contact-growth evaluation required.',
        accepted_time_steps=0,physical_providers=0,
        limits='Saved-state assessment only. New monitors do not rewrite prior accepted states. Growth ratios do not establish convergence order or a dynamical cause.')
    args.output.mkdir()
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:report[key] for key in ('status','all_three_times_available','available_density_passed','available_temporal_passed','automatic_extension_authorized','pending')}),flush=True)


if __name__=='__main__':main()

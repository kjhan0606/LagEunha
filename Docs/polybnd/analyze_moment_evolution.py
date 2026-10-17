"""Read saved evolution and independent geometric controls; no simulation."""
import json
import sys
from pathlib import Path
import numpy as np


def analyze(base,events,geometry):
    base,events,geometry=map(Path,(base,events,geometry))
    read=lambda p:[json.loads(s) for s in (p/'results.jsonl').read_text().splitlines()]
    first,extra,kin=read(base),read(events),read(geometry)
    assert len(first)==6 and len(extra)==6 and len(kin)==5
    assert all(r['status']=='complete' for r in first)
    smooth=sorted([r for r in first if r['n']==3],key=lambda r:-r['dt'])
    before=sorted([r for r in extra if r['end']==.006],key=lambda r:-r['dt'])
    rows=[r for r in first if r['n']==4]+[r for r in extra if r['end']==.02 and r['dt']<.00125]
    pairs=[]
    for g in kin:
        r=next(r for r in rows if r['dt']==g['dt'])
        pairs.append(dict(dt=g['dt'],fluid_pressure_error=r['final']['pressure_error'],
            kinematic_relative_area_drift=g['max_relative_area_drift'],predicted_pressure_error=.6*g['max_relative_area_drift']))
    def orders(rows):
        return [float(np.log2(a['final']['pressure_error']/b['final']['pressure_error'])) for a,b in zip(rows,rows[1:])]
    h=json.loads((events/'n4_dt0.00125_end0.02.history.json').read_text())
    transitions=[dict(time=r['time'],pressure_before=h[i-1]['pressure_error'] if i else 0.,pressure_after=r['pressure_error']) for i,r in enumerate(h) if r['topology_changed']]
    result=dict(smooth_pressure_orders=orders(smooth),before_event_pressure_orders=orders(before),
        geometry_comparison=pairs,topology_events=transitions,
        all_stage_realizability_pass=False,
        explanation='Negative residual covariance occurs in RK stages and some completed states. Execution completion is not a physical realizability pass.')
    (base/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    axes[0].semilogy([r['time'] for r in h],[r['pressure_error'] for r in h],marker='o')
    for r in transitions:axes[0].axvline(r['time'],color='tab:red',alpha=.4,ls='--')
    axes[0].set(title='4x4, dt=0.00125; dashed: neighbor change',xlabel='Time',ylabel='Max pressure error')
    axes[1].loglog([r['dt'] for r in pairs],[r['fluid_pressure_error'] for r in pairs],'o-',label='4x4 moment evolution')
    axes[1].loglog([r['dt'] for r in pairs],[r['predicted_pressure_error'] for r in pairs],'x--',label='0.6 x area drift, geometry only')
    axes[1].loglog([r['dt'] for r in smooth],[r['final']['pressure_error'] for r in smooth],'s-',label='3x3 moment evolution')
    axes[1].set(xlabel='dt',ylabel='Max pressure error',title='Final t=0.02');axes[1].legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.25)
    fig.savefig(base/'moment_geometry_diagnostic.png',dpi=180);fig.savefig(base/'moment_geometry_diagnostic.pdf');plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':analyze(*sys.argv[1:4])

"""Analyze saved shear histories; no evolution."""
import json
from pathlib import Path
import sys
import numpy as np


def analyze(directory):
    d=Path(directory); rows=[json.loads(s) for s in (d/'results.jsonl').read_text().splitlines()]
    output=[]
    for r in rows:
        if r['status']!='complete':output.append(r);continue
        hist=[json.loads(s) for s in (d/r['history']).read_text().splitlines()]
        initial=r['initial'];final=r['final']
        dp=np.diff([s['pressure_std'] for s in hist]);changes=np.array([s['added_edges']+s['removed_edges']>0 for s in hist[1:]])
        record=dict(n=r['n'],amplitude=r['amplitude'],closure=r['closure'],refine=r['refine'],steps=r['steps'],dt=r['dt'],
            initial_pressure_bias=r['initial_bias'],initial_pressure_std=initial['pressure_std'],
            final_pressure_std=final['pressure_std'],final_pressure_mean=final['pressure_mean'],
            final_projected_pressure_rms=final['pressure_rms_from_projected_exact'],
            final_transverse_velocity_rms=final['transverse_velocity_rms'],final_vx_error_rms=final['velocity_x_error_rms'],
            final_density_rms=final['density_rms'],mode_retention=r['mode_retention'],
            kinetic_change=r['kinetic_change'],internal_change=r['internal_change'],
            energy_relative=r['energy_relative'],mass_flux=r['mass_flux'],pressure_min=r['pressure_min'],area_min=r['area_min'],
            topology_steps=r['topology_steps'],edge_additions=r['edge_additions'],edge_removals=r['edge_removals'],
            max_pressure_std=max(s['pressure_std'] for s in hist),
            mean_abs_pressure_std_step_change_with_topology=float(np.mean(abs(dp[changes]))) if np.any(changes) else None,
            mean_abs_pressure_std_step_change_without_topology=float(np.mean(abs(dp[~changes]))) if np.any(~changes) else None,
            transverse_rms_history=[s['transverse_velocity_rms'] for s in hist],times=[s['time'] for s in hist])
        output.append(record)
    refinement=None
    candidates=[r for r in rows if r.get('status')=='complete' and r['n']==16 and r['amplitude']==.5 and r['closure']=='conservative']
    if len(candidates)==2:
        a,b=sorted(candidates,key=lambda r:r['refine']);ya=np.array(a['final_state']);yb=np.array(b['final_state'])
        refinement=dict(max_state_difference_by_column=np.max(abs(ya-yb),axis=0).tolist(),
            pressure_std_difference=abs(a['final']['pressure_std']-b['final']['pressure_std']),
            mode_retention_difference=abs(a['mode_retention']-b['mode_retention']))
    result=dict(cases=output,time_refinement=refinement)
    (d/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    for r in output:print({k:v for k,v in r.items() if k not in ('transverse_rms_history','times')})
    print('time_refinement',refinement)
    return result


if __name__=='__main__':analyze(sys.argv[1])


def plot_histories(directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d=Path(directory);rows=[json.loads(s) for s in (d/'results.jsonl').read_text().splitlines()]
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    for r in rows:
        if r['status']!='complete' or r['amplitude']!=.5 or r['refine']!=1:continue
        h=[json.loads(s) for s in (d/r['history']).read_text().splitlines()];t=[s['time'] for s in h]
        label=f'{r["n"]}x{r["n"]}, {r["closure"]}'
        color='tab:blue' if r['n']==8 else 'tab:orange';ls='-' if r['closure']=='conservative' else '--'
        axes[0,0].plot(t,[s['pressure_std']/.6*100 for s in h],ls,color=color,label=label)
        axes[0,1].plot(t,[s['transverse_velocity_rms']/.5 for s in h],ls,color=color)
        axes[1,0].plot(t,[s['mode_amplitude']/h[0]['mode_amplitude'] for s in h],ls,color=color)
        axes[1,1].plot(t,[(s['internal']-h[0]['internal'])/h[0]['kinetic'] for s in h],ls,color=color)
    titles=['Pressure spatial std / P0 (%)','Transverse velocity RMS / A','Shear mode / initial mode','Internal-energy gain / initial resolved KE']
    for ax,title in zip(axes.ravel(),titles):ax.set_title(title);ax.set_xlabel('Time');ax.grid(alpha=.25)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Fixed-mass ALE, sinusoidal shear A=0.5; short baseline test')
    fig.savefig(d/'shear_history.png',dpi=180);fig.savefig(d/'shear_history.pdf');plt.close(fig)

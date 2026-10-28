"""Plot completed local AV comparison; run on lageunha, not the login node."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap=argparse.ArgumentParser()
ap.add_argument('root',type=Path)
args=ap.parse_args()
data=json.loads((args.root/'analysis.json').read_text())
fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True)
summary={}
for name,case in data['cases'].items():
    h=case['history'];t=[r['time'] for r in h]
    gap=[r['pairs'][0]['gap_over_initial_h'] for r in h]
    speed=[r['regions']['inner']['circular_velocity_error_mass_RMS'] for r in h]
    radial=[r['regions']['inner']['radial_displacement_from_t0_mass_RMS'] for r in h]
    heat=[r['regions']['inner']['entropy_change_mass_mean'] for r in h]
    for ax,y in zip(axes.flat,[gap,speed,radial,heat]):ax.plot(t,y,label=name)
    k=min(range(len(gap)),key=gap.__getitem__)
    summary[name]={'sampled_min_gap_over_h':gap[k],'sampled_min_gap_time':t[k],
                   'final_time':t[-1],'final_gap_over_h':gap[-1],
                   'final_separation_rate':h[-1]['pairs'][0]['separation_rate'],
                   'final_inner':h[-1]['regions']['inner'],
                   'initial_inner':h[0]['regions']['inner'],
                   'accepted_steps':case['accepted_steps']}
for ax,label in zip(axes.flat,['803–804 separation / initial h','Inner circular velocity error (mass RMS)',
                              'Inner radial displacement from t0 (mass RMS)','Inner mean entropy change from restart']):
    ax.set_ylabel(label);ax.grid(alpha=.25);ax.legend()
for ax in axes[-1]:ax.set_xlabel('Absolute simulation time')
fig.suptitle('Local t6 restart: compression blend comparison (not a long-time validation)')
fig.tight_layout()
for ext in ['png','pdf']:
    target=args.root/f'comparison.{ext}'
    assert not target.exists(),target
    fig.savefig(target,dpi=160)
with (args.root/'comparison_summary.json').open('x') as f:json.dump(summary,f,indent=2);f.write('\n')
print(json.dumps(summary,indent=2))

"""Read explicit 2D batch paths only; summarize completed records and plot errors."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RUNS={
 'old':'Codex_B_adaptive_moments_20261010_000207',
 'compatible':'Codex_B_compatible_moments_20261010_000545',
 'baseline':'Codex_B_2d_state_suite_20261010_000859',
 'stress':'Codex_B_residual_stress_20261010_001101',
 'midpoint':'Codex_B_implicit_moments_20261010_001354',
 'linear_pressure':'Codex_B_linear_pressure_20261010_001601',
 'factor':'Codex_B_factor_moments_20261010_002319',
 'factor_2d':'Codex_B_factor_2d_20261010_002650',
 'kinetic':Path('/tmp/lagforce_kinetic_factors_path').read_text().strip().split('/')[-1]}
EXPECTED={'old':13,'compatible':13,'baseline':7,'stress':7,'midpoint':6,'linear_pressure':7,'factor':8,'factor_2d':6,'kinetic':8}
root=Path('/gpfs/kjhan/LagForce');data={};summary={}
for name,folder in RUNS.items():
    p=root/folder/'results.jsonl';rows=[]
    if p.exists():
        for line in p.read_text().splitlines():
            try:rows.append(json.loads(line))
            except json.JSONDecodeError:continue
    data[name]=rows;done=[r for r in rows if r['status']=='completed']
    summary[name]=dict(path=str(p.parent),records=len(rows),expected=EXPECTED[name],completed=len(done),failed=[r for r in rows if r['status']!='completed'],
        max_ledger_error=max((r['final']['ledger_error'] for r in done),default=None))
Path('bundle_summary.json').write_text(json.dumps(summary,indent=2))
fig,axes=plt.subplots(1,2,figsize=(10,4.2))
for name,label in [('baseline','Compatible pressure'),('stress','+ residual stress'),('linear_pressure','+ linear pressure')]:
    for ax,kind,key in [(axes[0],'wave','pressure_L1'),(axes[1],'gresho','velocity_L1')]:
        rows=[r for r in data[name] if r['status']=='completed' and r['kind']==kind and r.get('max_step')==.0025]
        rows=sorted(rows,key=lambda r:r['n'])
        if rows:ax.plot([r['n'] for r in rows],[r['errors'][key] for r in rows],'o-',label=label)
axes[0].set(title='Acoustic pressure error / amplitude',ylabel='Cell-average L1',xlabel='N (N x N cells)')
axes[1].set(title='Gresho velocity error (t=0.02)',ylabel='Cell-average L1',xlabel='N (N x N cells)')
for ax in axes:ax.grid(alpha=.3);ax.set_xticks([4,6,8]);ax.legend(fontsize=8)
fig.tight_layout();fig.savefig('convergence.png',dpi=180);fig.savefig('convergence.pdf');plt.close(fig)
print(json.dumps(summary,indent=2))

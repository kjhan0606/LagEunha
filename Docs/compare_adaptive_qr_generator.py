"""Compare saved adaptive one-step derivatives to the matching force generator."""
import json
from pathlib import Path
b=Path('/gpfs/kjhan/LagForce')
root=b/'Codex_adaptive_qr_steps_20261018'
g=json.loads((b/'Codex_adaptive_qr_periodic_probes_20261018/generator_analysis.json').read_text())
s=json.loads((root/'analysis.json').read_text())
assert abs(g['time']-s['time'])<1e-12 and s['adaptive_qr']
out={'time':s['time'],'regions':{},'scope':'Same saved states and adaptive force operator; forward finite difference versus instantaneous generator. Not a formal RK4 order, eigenmode, or long-time stability test. Tiny regional signals can be roundoff-sensitive.'}
for region,row in g['regions'].items():
 pred=row['total_generator'];comparisons={}
 for dt,values in s['comparisons'].items():
  value=values[region]['measured_forward_derivative']
  comparisons[dt]={'measured':value,'absolute_difference':value-pred,'relative_error':(value-pred)/abs(pred) if pred else None}
 out['regions'][region]={'generator':pred,'comparisons':comparisons}
with (root/'generator_comparison.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))

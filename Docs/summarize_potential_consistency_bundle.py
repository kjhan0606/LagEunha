"""Aggregate approved performance, polynomial and four-state diagnostics."""
import collections,hashlib,json
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');out=Path('Docs/canonical_potential_consistency_bundle_2026-10-18.json')
assert not out.exists(),'Preserve prior report'
paths=dict(performance=Path('Docs/canonical_moment_flux_provider_comparison_2026-10-18.json'),
    polynomial=base/'Codex_potential_polynomial_20261018/report.json',
    directional=base/'Codex_potential_directional_20261018/report.json')
reports={name:json.loads(p.read_text()) for name,p in paths.items()}
prior0=json.loads((base/'Codex_adjacency_prior_impact_20261018/midpoint/priors.json').read_text())['old']
transitions={}
for label in ('plus','minus','plus_half','minus_half'):
    priors=json.loads((base/'Codex_potential_directional_20261018'/(label+'_priors.json')).read_text())
    counts=collections.Counter();contrast=[]
    for p,q in zip(prior0,priors):
        old=(p['flat'],p.get('reason'),p.get('direction'),p['stencil_cells'])
        new=(q['flat'],q.get('reason'),q.get('direction'),q['stencil_cells'])
        if old!=new:counts[str((old,new))]+=1;contrast.append(p['high']/p['low']-1)
    transitions[label]=dict(counts=dict(counts),initial_contrast_minus_one=dict(minimum=min(contrast),median=float(np.median(contrast)),maximum=max(contrast)) if contrast else None)
result=dict(status='completed_bundle_for_Astra_evaluation',reports=reports,prior_branch_transitions=transitions,
    derivative_caution='Branch changes preclude claiming a confirmed global smooth derivative from this experiment; no additional epsilon runs were authorized or performed.',
    source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths.values()})
out.write_text(json.dumps(result,indent=2))
print(json.dumps(dict(full_step=reports['polynomial']['full_step'],derivatives=reports['directional']['derivatives'],branch_transitions=transitions),indent=2))

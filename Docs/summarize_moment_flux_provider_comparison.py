"""Aggregate the single approved scalar/batch accumulation comparison pair."""
import hashlib,json
from pathlib import Path
root=Path('/gpfs/kjhan/LagForce/Codex_moment_flux_provider_comparison_20261018')
output=Path('Docs/canonical_moment_flux_provider_comparison_2026-10-18.json')
assert not output.exists(),'Preserve prior result'
reports={name:json.loads((root/name/'report.json').read_text()) for name in ('scalar','batch')}
a=reports['scalar'];b=reports['batch']
assert a['status']==b['status']=='passed'
for key in ('seed_sha256','state_sha256','factors_sha256'):assert a['manifest'][key]==b['manifest'][key]
assert a['iterations']==b['iterations']==3
for r in reports.values():assert r['manifest']['flags']['LAGFORCE_BATCH_FACE_STATE']=='1'
result=dict(status='passed',reports=reports,
    speedup={key:a['timers']['provider_total'][key]/b['timers']['provider_total'][key] for key in ('wall_seconds','process_cpu_seconds')},
    paired_input_hashes={key:a['manifest'][key] for key in ('seed_sha256','state_sha256','factors_sha256')},
    source_sha256={str(root/name/'report.json'):hashlib.sha256((root/name/'report.json').read_bytes()).hexdigest() for name in reports},
    limits='One provider evaluation each, face-state batch on for both. Shared-host observation, no repeated tuning, no evolution. The copied per-report profile_based_ideal_Amdahl_limit belongs to the earlier density/velocity operator and is NOT applicable to accumulation; ignore that inherited field.')
output.write_text(json.dumps(result,indent=2))
print(json.dumps(dict(status=result['status'],speedup=result['speedup'],batch_comparisons=b['comparisons']),indent=2))

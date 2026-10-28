"""Single-pair face-record performance and numerical-gate aggregation."""
import hashlib,json
from pathlib import Path
root=Path('/gpfs/kjhan/LagForce/Codex_face_record_provider_comparison_20261019')
out=Path('Docs/canonical_face_record_provider_comparison_2026-10-19.json')
assert not out.exists(),'Preserve prior comparison'
reports={name:json.loads((root/name/'report.json').read_text()) for name in ('scalar','batch')}
a=reports['scalar'];b=reports['batch'];assert a['status']==b['status']=='passed'
for key in ('seed_sha256','state_sha256','factors_sha256'):assert a['manifest'][key]==b['manifest'][key]
for r in reports.values():
 assert r['iterations']==3
 for key in ('LAGFORCE_BATCH_FACE_STATE','LAGFORCE_VECTOR_MOMENT_FLUX'):assert r['manifest']['flags'][key]=='1'
record_totals={name:{metric:sum(row['record_total'][metric] for row in r['rhs_calls']) for metric in ('wall_seconds','process_cpu_seconds')} for name,r in reports.items()}
result=dict(status='passed',reports=reports,record_totals=record_totals,
 speedup={metric:a['timers']['provider_total'][metric]/b['timers']['provider_total'][metric] for metric in ('wall_seconds','process_cpu_seconds')},
 export_excluded_speedup={metric:a['excluding_face_export'][metric]/b['excluding_face_export'][metric] for metric in ('wall_seconds','process_cpu_seconds')},
 source_sha256={str(root/name/'report.json'):hashlib.sha256((root/name/'report.json').read_bytes()).hexdigest() for name in reports},
 limits='One call each, shared-host observation. Scalar record total includes exact-star calls inside original scalar loop; native subregions separate preparation/star/packing. Export-excluded timing subtracts only the face export, not all diagnostic work. No closure change or new evolution in this comparison.')
out.write_text(json.dumps(result,indent=2))
print(json.dumps({k:result[k] for k in ('status','speedup','export_excluded_speedup','record_totals')},indent=2))

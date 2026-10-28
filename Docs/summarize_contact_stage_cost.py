"""Measured existing stage wall times; no timing experiment or rerun."""
from pathlib import Path
import json,statistics
from datetime import datetime,timezone
base=Path('/gpfs/kjhan/LagForce');reports=[]
for name,dt in [('Codex_contact_prior_evolution_dt1e4_20261018',1e-4),('Codex_contact_prior_evolution_dt5e5_20261018',5e-5)]:
 root=base/name;path=root/'events.jsonl'
 with path.open('rb') as stream:
  stream.seek(0,2);size=stream.tell();stream.seek(max(0,size-1000000));data=stream.read()
 events=[]
 for line in data.splitlines():
  try:events.append(json.loads(line))
  except (ValueError,UnicodeDecodeError):pass
 stages=[r for r in events if 'stage_seconds' in r];seconds=[r['stage_seconds'] for r in stages];accepted=[r['accepted'] for r in events if 'accepted' in r];last=accepted[-1] if accepted else None;summary=root/'summary.json'
 median=statistics.median(seconds);step=2*median
 reports.append(dict(root=str(root),dt=dt,latest_accepted=last,completed_summary=summary.exists(),recorded_stages=stages,measured_stage_count=len(stages),stage_wall_seconds=dict(minimum=min(seconds),median=median,maximum=max(seconds),sum=sum(seconds)),cost_projection_patch_only=[dict(horizon=t,assumed_dt=dt,steps=round(t/dt),assumed_stages_per_step=2,wall_seconds=step*round(t/dt),wall_days=step*round(t/dt)/86400) for t in (1.,10.)],limits='Observed wall times under concurrent server load, not exclusive CPU time or isolated per-component profiling. Cost projection only at same1590-cell patch and fixed dt; no full-domain scaling guarantee. Two stages per midpoint step assumed. No timestep increase approved.'))
report=dict(observed_utc=datetime.now(timezone.utc).isoformat(),scope=__doc__,runs=reports,component_timing='Existing log has total provider stage time only. Geometry/density/reconstruction/root/face/gravity bottlenecks cannot be attributed individually from these fields.')
print(json.dumps(report,indent=2))

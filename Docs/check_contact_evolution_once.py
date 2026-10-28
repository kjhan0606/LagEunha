"""One bounded progress read; pair with a live session/process observation.

Run at most once per minute. A summary or log alone is not a live-process test.
Never launches, stops or restarts a calculation.
"""
from datetime import datetime,timezone
import json
from pathlib import Path
base=Path('/gpfs/kjhan/LagForce');rows=[]
for name in ('Codex_contact_prior_evolution_dt1e4_20261018','Codex_contact_prior_evolution_dt5e5_20261018'):
    root=base/name
    with (root/'events.jsonl').open('rb') as stream:
        stream.seek(0,2);size=stream.tell();stream.seek(max(0,size-100000));data=stream.read()
    events=[]
    for line in data.splitlines():
        try:events.append(json.loads(line))
        except (json.JSONDecodeError,UnicodeDecodeError):pass
    accepted=[r['accepted'] for r in events if 'accepted' in r]
    summary=root/'summary.json'
    rows.append(dict(root=str(root),latest_accepted=accepted[-1] if accepted else None,
        latest_stage_event=events[-1] if events else None,
        summary=json.loads(summary.read_text()) if summary.exists() else None))
print(json.dumps(dict(observed_utc=datetime.now(timezone.utc).isoformat(),runs=rows,scope=__doc__),indent=2))

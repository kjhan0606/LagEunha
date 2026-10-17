"""Event-driven final scientific report after adaptive computation exits."""
import json
from pathlib import Path
from datetime import datetime,timezone
RUN=Path('/gpfs/kjhan/LagForce/Codex_P_adaptive_coupled_H200_20261010_200459')
OUTPUT=Path('/gpfs/kjhan/LagForce/Codex_P_monitor_20261010_214343')
report=dict(checked_at=datetime.now(timezone.utc).isoformat(),cases={})
for name in ('shear','wave'):
    p=RUN/name/'result.json'
    if not p.exists():
        report['cases'][name]=dict(status='incomplete',passed=False)
        continue
    r=json.loads(p.read_text())
    passed=r.get('status')=='completed' and all(r.get(k) is True for k in
        ('state_gate','mass_path_gate','generator_path_gate'))
    report['cases'][name]=dict(result=r,passed=passed)
report['passed']=all(r['passed'] for r in report['cases'].values())
report['next_action']='review coupling results and select next pre-Kepler validation' if report['passed'] else 'inspect failure or resume incomplete checkpoint; do not declare success'
(OUTPUT/'completion_report.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report))
raise SystemExit(0 if report['passed'] else 2)

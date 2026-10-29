"""Audit both same-GPU impact cases and full native-record parity."""
import json
from pathlib import Path
import sys
import numpy as np
from audit_kepler_wall_events_20261029 import analyze
from audit_kepler_restart_bytes_20261029 import compare
from kepler_native_global_ids_20261029 import ordered_global

root = Path(sys.argv[1])
ids = 2*np.arange(55351, dtype=np.uint64)+7
results = {}
for name in ('impact_A', 'impact_B'):
    results[name] = analyze(root/name, reader=lambda p: ordered_global(p, ids))
endpoints = [Path(results[name]['endpoint']) for name in ('impact_A', 'impact_B')]
parity = compare(*endpoints)
checks = {'both_physical_gates': all(r['passed'] for r in results.values()),
          'both_real_impacts': all(len(r['impacts']) > 0 for r in results.values()),
          'both_real_rollbacks': all(r['rollbacks'] > 0 for r in results.values()),
          'full_record_footer_parity': parity['passed']}
result = {'passed': all(checks.values()), 'checks': checks, 'cases': results, 'parity': parity,
          'scope': 'One real late Q2 collision with rollback in each native candidate; no long-time, MPI or convergence claim.'}
with (root/'impact_audit.json').open('x') as stream:
    json.dump(result, stream, indent=2); stream.write('\n')
print(json.dumps(result, indent=2))
if not result['passed']: raise SystemExit(1)

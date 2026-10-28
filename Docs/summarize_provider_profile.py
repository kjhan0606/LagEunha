"""Summarize one completed provider profile; no new replay or tuning."""
import argparse
import hashlib
import json
import math
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--report', type=Path, required=True)
args = parser.parse_args()
d = json.loads(args.report.read_text())
assert d['status'] == 'completed'
assert d['iterations'] == d['original_iterations']
assert all(row['maximum_tolerance_fraction'] <= 1 for row in d['comparisons'].values())
total = d['timers']['provider.total']
wall = sum(row['exclusive_wall_seconds'] for row in d['timers'].values())
cpu = sum(row['exclusive_cpu_seconds'] for row in d['timers'].values())
assert math.isclose(wall, total['wall_seconds'], rel_tol=1e-12, abs_tol=1e-8)
assert math.isclose(cpu, total['cpu_seconds'], rel_tol=1e-12, abs_tol=1e-8)
rows = []
for name, row in d['timers'].items():
    rows.append(dict(name=name, **row,
        fraction_provider_wall_exclusive=row['exclusive_wall_seconds']/total['wall_seconds'],
        fraction_provider_cpu_exclusive=row['exclusive_cpu_seconds']/total['cpu_seconds']))
rows.sort(key=lambda row: row['exclusive_cpu_seconds'], reverse=True)
print(json.dumps(dict(scope=__doc__, report=str(args.report),
    report_sha256=hashlib.sha256(args.report.read_bytes()).hexdigest(),
    numerical_reproduction=d['comparisons'], outer_iterations=d['iterations'],
    provider_total=total, exclusive_sum_wall_seconds=wall, exclusive_sum_cpu_seconds=cpu,
    components_sorted_by_exclusive_process_cpu=rows,
    concurrency=d['concurrency'], environment=d['environment'],
    limits='One instrumented replay with concurrent fine load and timer overhead. Exclusive clock spans remove measured children, not uninstrumented work. Ranking is measurement evidence for Astra, not an optimization choice, speedup result or long-term physics pass.'), indent=2))

"""Summarize existing per-cell preparation records without evaluating physics."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cells', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    digest = hashlib.sha256()
    groups = {'all': [], 'flat_prior': [], 'nonflat_prior': []}
    with args.cells.open('rb') as stream:
        for line in stream:
            digest.update(line)
            row = json.loads(line)
            assert row['status'] == 'passed'
            groups['all'].append(row)
            groups['flat_prior' if row['prior']['flat'] else 'nonflat_prior'].append(row)
    summaries = {}
    for name, rows in groups.items():
        if not rows:
            summaries[name] = {'cells': 0}
            continue
        item = {'cells': len(rows)}
        for key in ('nodes24', 'nodes48'):
            values = np.array([row[key] for row in rows], dtype=np.int64)
            index = int(values.argmax())
            item[key] = dict(total=int(values.sum()), median=float(np.median(values)),
                             p95=float(np.quantile(values, .95)), maximum=int(values[index]),
                             maximum_id=rows[index]['id'])
        keys = sorted(set.intersection(*(set(row['timers']) for row in rows)))
        item['summed_cell_wall_seconds'] = {
            key: float(sum(row['timers'][key] for row in rows)) for key in keys}
        assert all(np.isfinite(value) and value >= 0
                   for value in item['summed_cell_wall_seconds'].values())
        summaries[name] = item
    report = dict(status='completed_saved_density_preparation_cost_summary',
                  source=str(args.cells.resolve()), source_sha256=digest.hexdigest(),
                  groups=summaries, physical_providers=0, accepted_time_steps=0,
                  limits=['Summed per-cell wall times overlap across workers; these are not job wall times.',
                          'Independent q48 is checkpoint validation, not every-RHS evolution cost.',
                          'Node counts alone do not justify quadrature-order reduction.'])
    args.output.mkdir()
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()

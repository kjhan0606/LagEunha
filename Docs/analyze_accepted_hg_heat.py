"""Audit HG heat using accepted steps only (canonical MPI1 split runs).

Regional HG bins use instantaneous radius, not material membership. The
remaining global internal-energy change includes reversible work and all
non-HG updates; it must not be labelled TCF dissipation.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def analyze(raw):
    pending = []
    initial = previous = None
    records = []
    rejected = 0
    max_mismatch = 0.0
    def fields(line):
        return dict(token.split('=', 1) for token in line.split()[1:] if '=' in token)
    for line in raw.decode().splitlines():
        if line.startswith('[HG_TRANSFER]'):
            d = fields(line)
            assert int(d['rank']) == 0, 'Only MPI1 audit supported'
            pending.append(d)
        elif line.startswith('[RETRY] step rejected'):
            pending.clear()
            rejected += 1
        elif line.startswith('[ACCEPTED_INITIAL]'):
            initial = previous = float(fields(line)['U'])
        elif line.startswith('[ACCEPTED_STATE]'):
            d = fields(line)
            assert previous is not None
            assert int(d['step']) == len(records) + 1
            assert [int(v['bin']) for v in pending] == [0, 1, 2, 3] * 2
            assert all(abs(2*float(v['tau'])-float(d['dt'])) <= 1e-9*float(d['dt']) for v in pending)
            heat = [sum(float(v['dU']) for v in pending if int(v['bin']) == b) for b in range(4)]
            mismatch = abs(sum(heat)-float(d['hg_dU']))
            assert math.isfinite(mismatch) and mismatch < 1e-11
            max_mismatch = max(max_mismatch, mismatch)
            du = float(d['U'])-previous
            records.append({'step':int(d['step']), 't':float(d['elapsed']),
                            'delta_U':du, 'HG_heat_by_bin':heat,
                            'non_HG_delta_U':du-sum(heat)})
            previous = float(d['U'])
            pending.clear()
    assert records
    windows = []
    for lo, hi in [(0,2),(2,4),(4,6)]:
        selected = [r for r in records if lo < r['t'] <= hi]
        if not selected:
            continue
        windows.append({'endpoint_selection': [lo,hi], 'last_t':selected[-1]['t'],
                        'steps':len(selected),
                        'delta_U':sum(r['delta_U'] for r in selected),
                        'HG_heat_by_bin':[sum(r['HG_heat_by_bin'][b] for r in selected) for b in range(4)],
                        'non_HG_delta_U':sum(r['non_HG_delta_U'] for r in selected)})
    return {'accepted_steps':len(records), 'last_t':records[-1]['t'],
            'rejected_trials_discarded':rejected,
            'uncommitted_HG_records_ignored':len(pending),
            'max_HG_ledger_mismatch':max_mismatch, 'windows':windows,
            'scope':'MPI1; instantaneous radial bins; windows select step endpoints, not interpolated exact times. Non-HG delta U includes compression and cannot be interpreted as dissipative heating.'}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('log')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    # Bound the read of a running log to the current file size and complete lines.
    with open(args.log, 'rb') as f:
        size = Path(args.log).stat().st_size
        raw = f.read(size)
    raw = raw[:raw.rfind(b'\n')+1]
    result = analyze(raw)
    result.update(log=args.log, prefix_bytes=len(raw), prefix_sha256=hashlib.sha256(raw).hexdigest())
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))

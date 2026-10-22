#!/usr/bin/env python3
"""Read-only integrity audit. A valid rejection is never a physics pass."""
import csv
import hashlib
import json
from pathlib import Path
import struct
import sys

ROOT = Path(__file__).resolve().parent


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(data)
    return digest.hexdigest()


def main():
    run = Path(sys.argv[1]).resolve()
    d = json.loads((run / 'contract.json').read_text())
    preflight = json.loads((run / 'duffy_preflight/contract.json').read_text())
    assert preflight['duffy_AD2_contract_pass']
    assert preflight['sources'] == 384
    assert preflight['enclosure_mismatches'] == 0
    assert preflight['scalar_comparisons'] == 295296
    assert preflight['min_positive_J'] > 0
    assert preflight['whole_real_other_denominator_lower'] > 0
    contributions = rows(run / 'source_target_contributions.csv')
    assert len(contributions) == 384
    assert {int(r['source']) for r in contributions} == set(range(384))
    counts = [sum(int(r['classification']) == k for r in contributions)
              for k in range(3)]
    assert counts == [344, 4, 36]
    failed = sum(r['all12_complete'] != '1' for r in contributions)
    geom_failed = sum(r['geometry_complete'] != '1' for r in contributions)
    assert failed == d['failed_all12_sources']
    assert geom_failed == d['unproved_actual_source_domains']
    assert d['CUDA_reader_original_points'] == 817152
    budgets = rows(run / 'all12_whole_target_budget.csv')
    assert len(budgets) == 12
    assert len({r['moment'] for r in budgets}) == 12
    if failed:
        assert not d['target_published']
        assert not d['all12_whole_target_budgets_pass']
        assert all(r['total_error'] == 'unproved' for r in budgets)
        assert all(r['budget_pass'] == '0' for r in budgets)
        assert d['caller_target_unchanged_on_rejection']
        packet = (run / 'caller_target_packet.bin').read_bytes()
        assert struct.unpack_from('=QQ', packet) == (9, 1234)
        assert (run / 'native_exit_code.txt').read_text().strip() == '2'
    original = ROOT / 'shared_trace_rhs_v1_mpi4/shared_state.bin'
    expected = sha(original)
    assert expected == '74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005'
    for relative in ('caller_unchanged_state.bin',
                     'immutable_heat_generating_state.bin',
                     'duffy_preflight/caller_unchanged_state.bin'):
        assert sha(run / relative) == expected
    context = rows(run / 'frozen_geometric_EOS_context.csv')
    assert len(context) == 24
    assert len({r['cell'] for r in context}) == 24
    reasons = {}
    for r in contributions:
        if r['all12_complete'] != '1':
            reasons[r['reason']] = reasons.get(r['reason'], 0) + 1
    output = {
        'scope': 'read-only audit of one native whole-target call',
        'input_run': str(run),
        'integrity_audit_pass': True,
        'actual_native_result': d,
        'duffy_value_first_second_derivative_contract': preflight,
        'all384_classification': counts,
        'original_state_and_immutable_heat_generator_SHA256': expected,
        'failed_source_reasons': reasons,
        'old_problematic_sources': {
            r['source']: r for r in contributions
            if r['source'] in ('118', '246', '254')
        },
        'whole_all12_physical_certificate_pass': d['target_published'],
        'long_Kepler_goal_achieved': False,
    }
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()

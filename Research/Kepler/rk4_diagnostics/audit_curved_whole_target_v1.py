"""Exact-decoded ONE finite-packet target failure ledger; no state changes."""
import csv
import hashlib
import json
from collections import Counter
from decimal import localcontext
from fractions import Fraction as F
from pathlib import Path
from audit_curved_inverse_physical_v1 import rbin, r80, rounded80, string

ROOT = Path(__file__).resolve().parent


def read(folder, name):
    return list(csv.DictReader((folder / name).open()))


def main():
    with localcontext() as ctx:
        ctx.prec = 85
        four = ROOT / 'curved_whole_target_v3_mpi4'
        one = ROOT / 'curved_whole_target_v3_mpi1'
        contracts = [json.loads((p / 'native_contract.json').read_text()) for p in (four, one)]
        for p, contract in zip((four, one), contracts):
            assert (p / 'native_exit_code.txt').read_text().strip() == '2'
            assert contract['all384_source_rows'] == 384
            assert contract['geometry_complete'] and contract['source_pairs_examined'] == 73536
            assert not contract['sources_complete'] and contract['failed_sources'] == 3
            assert not contract['published_one_target']
            assert contract['epoch_reject_preserves_output']
            assert contract['missing_row_reject_preserves_output']
            assert contract['actual_failure_preserves_output']
            assert not (p / 'validated_target.bin').exists()
        same_files = ('source_contributions.csv', 'whole_target_ledger.csv',
                      'canonical_seams.csv', 'canonical_vs_duplicated_seams.csv',
                      'nonadjacent_pairs.csv', 'external_edges.csv', 'caller_after.bin')
        reproducible = {name: (four / name).read_bytes() == (one / name).read_bytes()
                        for name in same_files}
        assert all(reproducible.values())
        rows = read(four, 'source_contributions.csv')
        assert len(rows) == 384 * 12
        by = {(int(r['source']), r['moment']): r for r in rows}
        mass_rows = [r for r in rows if r['moment'] == 'mass']
        failures = [{k: r[k] for k in ('source', 'cell', 'classification', 'reason')}
                    for r in mass_rows if r['complete'] == '0']
        assert [int(r['source']) for r in failures] == [118, 246, 254]
        lower = next(r for r in read(ROOT / 'curved_target_lower_v1_225', 'target_lower.csv')
                     if r['target'] == '68526')
        U = r80(lower['U_lower'])
        ledger = read(four, 'whole_target_ledger.csv')
        result = []
        for row in ledger:
            moment = row['moment']
            # This is a PARTIAL candidate and nonnegative known-bound sum.
            # Unknown errors of three sources remain unknown, never zero.
            partial, rounding, analytic, field = F(0), F(0), F(0), F(0)
            for s in range(384):
                entry = by[(s, moment)]
                if entry['complete'] != '1':
                    assert entry['native'] == entry['analytic_remainder'] == 'unproved'
                    continue
                value = r80(entry['native'])
                exact = partial + value
                rounded = rounded80(exact)
                rounding += abs(rounded - exact)
                partial = rounded
                analytic += rbin(entry['analytic_remainder'])
                field += r80(entry['native_geometry_nodes_sum_error'])
            assert partial == r80(row['candidate'])
            assert rounding <= rbin(row['final_sum_rounding'])
            assert field <= rbin(row['native_geometry_node_source_sum_error'])
            seam = rbin(row['canonical_seam_bound'])
            cap = rbin(row['vertex_cap_bound'])
            tolerance = rbin(row['unchanged_physical_tolerance'])
            if moment in ('K', 'U', 'potential', 'thermal_rate'):
                exact_tau = U / (10**9 if moment == 'thermal_rate' else 10**10)
                assert 0 < tolerance <= exact_tau
                assert exact_tau - tolerance < exact_tau * F(1, 10**32)
            known = field + analytic + seam + cap + rbin(row['final_sum_rounding'])
            assert row['total_absolute_error'] == row['error_over_tolerance'] == 'unproved'
            result.append({'moment': moment, 'partial_candidate': string(partial),
                           'completed_source_field_node_geometry_errors': string(field),
                           'completed_source_integral_errors': string(analytic),
                           'canonical_once_seam_bound': string(seam),
                           'vertex_cap_bound': string(cap),
                           'actual_binary80_sum_rounding': string(rounding),
                           'fixed_whole_target_tolerance': string(tolerance),
                           'sum_of_KNOWN_bound_terms_only': string(known),
                           'known_bound_sum_over_tolerance': string(known / tolerance),
                           'seam_bound_alone_over_tolerance': string(seam / tolerance),
                           'full_absolute_error': 'UNPROVED: three missing source certificates',
                           'pass': False})
        potential = next(r for r in result if r['moment'] == 'potential')
        assert F(potential['seam_bound_alone_over_tolerance']) > 1
        source_old = ROOT / 'curved_whole_target_v2_mpi4' / 'source_contributions.csv'
        assert source_old.read_bytes() == (four / 'source_contributions.csv').read_bytes()
        report = {
            'scope': 'fixed original384-source /24-cell finite packet target68526; rejected native call',
            'reused_original_source_certificates_without_reintegrating': True,
            'original_source_certificate_sha256': hashlib.sha256(source_old.read_bytes()).hexdigest(),
            'classifications_including_failed': dict(Counter(r['classification'] for r in mass_rows)),
            'completed_sources': 381, 'certified_positive_points': 81920,
            'failed_source_certificates': failures,
            'all73536_source_pairs_examined': True,
            'relevant_pair_certificates': dict(Counter(r['certificate'] for r in read(four, 'nonadjacent_pairs.csv'))),
            'canonical_seam_parts': contracts[0]['canonical_target_seam_parts'],
            'external_reference_edges': 80, 'conservative_possible_target_external_edges': 4,
            'external_policy': 'exact finite source pushforward boundary; no exterior fill or global disk claim',
            'MPI1_MPI4_bitwise_reproducible_files': reproducible,
            'native_failure_preserves_caller_output': True,
            'old_equal_work_share_failures_preserved': True,
            'physical_U_lower_or_energy_tolerances_changed': False,
            'exact_global_C0_map_created': False,
            'actual_physical_error_over_tolerance_proved': False,
            'validated_target_published': False,
            'CUDA_run_or_physical_time_evolution': False,
            'full_remap_or_long_Kepler_stability_proved': False,
            'stop_old_broken_packet_promotion': True,
            'moments': result,
        }
        with (ROOT / 'CURVED_WHOLE_TARGET_LEDGER_20261021.json').open('x') as f:
            json.dump(report, f, indent=2)
        print('381/384 complete; target rejected; MPI1/4 retained-proof replay matches bitwise.')
        print('potential seam/tau', potential['seam_bound_alone_over_tolerance'])
        print('potential KNOWN bound terms/tau', potential['known_bound_sum_over_tolerance'])


if __name__ == '__main__':
    main()

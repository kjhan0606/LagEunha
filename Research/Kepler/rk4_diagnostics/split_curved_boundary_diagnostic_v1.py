"""Read-only decomposition of the retained boundary witness.

The recovered per-panel field envelope is a diagnostic quotient of rounded
CSV values, NOT a new certified upper bound or replacement gate.
"""
import csv
import json
from collections import defaultdict
from decimal import Decimal as D, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read(path):
    with (ROOT / path).open() as handle:
        return list(csv.DictReader(handle))


def main():
    summary = json.loads((ROOT / 'CURVED_ORIGINAL_ERROR_LEDGER_20261021.json').read_text())
    report = {'scope': 'diagnostic attribution of saved CSV bounds, not a new proof', 'cases': {}}
    with localcontext() as ctx:
        ctx.prec = 70
        for case in ('225', '137', 'lens'):
            worst = summary['cases'][case]['worst_partial_row']
            target, moment = worst['target'], worst['moment']
            sites = [101, 102, 103, 104] if case == 'lens' else [int(r['id']) for r in read('curved_map_algebra_v2/advected_sites.csv')]
            j = sites.index(target)
            source = int(read('curved_fold_transfer_v5_' + case + '/target_contribution_bounds.csv')[0]['source'])
            # The three present worst rows are all incoming receivers.
            assert source != target
            receivers = {(int(r['edge']), int(r['target'])): bool(int(r['possible']))
                         for r in read('curved_fold_transfer_v5_' + case + '/source_boundary_receivers.csv')}
            geometry = read('curved_error_ledger_v3_' + case + '/regular_boundary_forward.csv')
            selected, old_area = [], defaultdict(D)
            for row in geometry:
                edge = int(row['source_edge'])
                relevant = receivers[(edge, target)] if edge >= 0 else j in (int(row['i']), int(row['j']))
                if relevant:
                    selected.append(row)
                    old_area[int(row['panel'])] += D(row['old_reference_measure'])
            old_field = {}
            for row in read('curved_fold_transfer_v5_' + case + '/panel_moments.csv'):
                panel = int(row['panel'])
                if int(row['target']) == j and row['moment'] == moment and old_area[panel] > 0:
                    old_field[panel] = D(row['analytic_bound']) / old_area[panel]
            split = defaultdict(D)
            records = []
            for row in selected:
                panel = int(row['panel'])
                if panel not in old_field:
                    continue
                contribution = D(row['new_reference_measure']) * old_field[panel]
                method = int(row['method'])
                if method in (2, 3, 4):
                    delta = D(row['new_delta'])
                    for col in ('actual_graph_error', 'conversion_error', 'root_guard'):
                        split[col] += contribution * D(row[col]) / delta
                else:
                    split['preserved_axis_or_source_fold'] += contribution
                records.append({**row, 'diagnostic_moment_contribution': str(contribution)})
            other = sum((D(v) for k, v in worst['known_error_model_components'].items() if k != 'regular_boundary'), D(0))
            remaining = D(worst['fixed_incidence_partial_allocation']) - other
            boundary = D(worst['known_error_model_components']['regular_boundary'])
            required_ratio = max(D(0), remaining) / boundary
            dominant = max(records, key=lambda r: D(r['diagnostic_moment_contribution']))
            report['cases'][case] = {
                'target': target, 'moment': moment,
                'diagnostic_split': {k: str(v) for k, v in split.items()},
                'dominant_saved_boundary_row': dominant,
                'required_uniform_boundary_factor_at_most_ignoring_missing_source_gluing_and_quad': str(required_ratio),
                'dominant_delta_if_every_boundary_scales_uniformly_at_most': str(D(dominant['new_delta']) * required_ratio),
                'required_factor_is_a_design_target_not_an_error_bound': True,
            }
            print(case, 'split', {k: str(v) for k, v in split.items()},
                  'required-factor<=', required_ratio, 'dominant-panel', dominant['panel'])
    with (ROOT / 'CURVED_BOUNDARY_DIAGNOSTIC_20261021.json').open('x') as handle:
        json.dump(report, handle, indent=2)


if __name__ == '__main__':
    main()

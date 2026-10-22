"""Adopt the unchanged slot3203 plan after independent coverage and Astra review."""
from collections import Counter, defaultdict
from pathlib import Path
import json
import struct
from adopt_shared_ray_plans_v1 import write_csv
from audit_shared_regular_q_transfer_v1 import rows, sha
ROOT = Path(__file__).resolve().parent

def main():
    previous = ROOT/'shared_atlas_adoption_v4'
    coverage = ROOT/'shared_weighted_coverage_v1'
    run = ROOT/'shared_weighted_ray_v1_mpi1'
    audit = ROOT/'SHARED_WEIGHTED_RAY_20261021.json'
    review = ROOT/'ASTRA_SHARED_WEIGHTED_RAY_END_REVIEW_20261021.md'
    dest = ROOT/'shared_atlas_adoption_v5'
    assert not dest.exists()
    for directory in (previous, coverage, run):
        for name in ('input.sha256', 'output.sha256'):
            for line in (directory/name).read_text().splitlines():
                digest, path = line.split('  ', 1)
                assert sha(Path(path)) == digest
    assert 'slot3203' in review.read_text() and '채택할 수 있다' in review.read_text()
    proof = json.loads((coverage/'contract.json').read_text())
    assert proof['domain_coverage_and_native_geometry_density_clearance_complete']
    assert json.loads(audit.read_text())['independent_weighted_ray_accounting_pass']
    [plan] = proof['adoptable_slot_plans']
    assert plan['slot_id'] == 3203 and plan['panels'] == 304
    packet = run/'unchanged_original_plan_rank0.bin'
    data = packet.read_bytes(); header = struct.unpack_from('<8Q', data)
    assert header[3:5] == (704, 304) and len(data) == 64+304*704
    for j in range(304):
        assert struct.unpack_from('<4Q', data, 64+704*j) == (3203,133,j,68000)
    ledger = rows(previous/'effective_slot_plan_references.csv')
    costs = defaultdict(Counter); pending = []
    for row in ledger:
        slot = int(row['slot_id'])
        row['required_angular_proof_evaluator'] = 'original_angular'
        row['radius_branch_proof_reference'] = ''
        if slot == 3203:
            assert row['atlas_complete'] == '0' and row['active_panels'] == '0'
            row.update(atlas_complete='1', active_panels='304', complete_planned_queries='77824',
                provenance='adopted_unchanged_actual_plan_weighted_ray_germ_clearance',
                typed_packet=str(packet), typed_record_first='0', typed_record_count='304',
                required_reference_radius_evaluator='collected_origin_quadratic_positive_real_anchored_germ',
                required_angular_proof_evaluator='weighted_angular_same_exact_function',
                radius_branch_proof_reference=str(run/'new_weighted_axis_proofs_rank0.bin'))
        value = costs[int(row['target'])]
        value['pending'] += row['atlas_complete'] == '0'
        value['retained' if row['provenance']=='retained_original_target_plan' else 'new'] += int(row['active_panels'])
        if row['atlas_complete']=='0':
            pending.append(dict(slot_id=slot,source=int(row['source']),target=int(row['target']),original_failure=row['original_reason']))
    policy = ROOT/'SHARED_ATLAS_BUDGET_POLICY_20261021.json'
    cap = json.loads(policy.read_text())['new_per_target_query_limit']
    targets = []
    for target,value in sorted(costs.items()):
        known = (value['new']+value['retained'])*256
        targets.append(dict(target=target,incomplete_atlas_slots=value['pending'],active_new_panels=value['new'],
            retained_original_target_panels=value['retained'],known_positive_query_lower_bound=known,
            complete_planned_queries='UNKNOWN' if value['pending'] else known,
            authorized_query_cap=cap,complete_cost_pass=int(not value['pending'] and known<=cap),analytic_certificates_complete=0))
    assert len(ledger)==9216 and len(pending)==29 and sum(t['active_new_panels'] for t in targets)==34955
    assert next(t for t in targets if t['target']==68000)['known_positive_query_lower_bound']==580864
    helpers = [ROOT/name for name in ('material_shared_collected_radius.h','material_shared_weighted_angular.h','material_shared_origin_quadratic.h')]
    inputs = [Path(__file__).resolve(),previous/'contract.json',previous/'effective_slot_plan_references.csv',coverage/'contract.json',audit,review,packet,run/'new_weighted_axis_proofs_rank0.bin',policy]+helpers
    dest.mkdir()
    write_csv(dest/'effective_slot_plan_references.csv',ledger)
    write_csv(dest/'effective_target_cost.csv',targets)
    write_csv(dest/'remaining29_atlas_slots.csv',pending)
    result = dict(scope='slot3203 actual unchanged atlas and geometry/density clearance adopted after independent audit and completed Astra branch review',
        newly_completed_atlas_slots=[3203],newly_adopted_panels=304,active_new_panels=34955,
        retained_original_target_panels=1047,remaining_incomplete_atlas_slots=29,authorized_target_query_limit=cap,
        historical_pending_review_artifact_preserved=True,completed_branch_review=str(review),
        required_proof_helpers={str(p):sha(p) for p in helpers},
        real_branch_anchor='positive sqrt of P(t) at each saved real panel center; r=s*sqrt(P), no conjugation',
        future_analytic_evaluator_must_consume_required_helpers=True,
        new_physical_quadrature_queries=0,all_target_actual_atlas_complete=False,
        all_target_analytic_certificates_complete=False,full_quadrature_allowed=False,
        target_cost=targets,remaining_atlas_slots=pending)
    out=dest/'contract.json';out.write_text(json.dumps(result,indent=2)+'\n')
    (dest/'input.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in inputs))
    (dest/'output.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in (out,dest/'effective_slot_plan_references.csv',dest/'effective_target_cost.csv',dest/'remaining29_atlas_slots.csv')))
    print(json.dumps({k:result[k] for k in ('newly_completed_atlas_slots','active_new_panels','remaining_incomplete_atlas_slots','full_quadrature_allowed')},indent=2))

if __name__=='__main__': main()

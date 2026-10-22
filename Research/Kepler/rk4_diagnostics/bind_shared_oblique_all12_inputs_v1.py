"""Bind existing complex field snapshots; no native query or quadrature."""
from collections import Counter
from pathlib import Path
import csv
import hashlib
import json
import struct
from audit_shared_atlas_witness_v1 import wire
from audit_shared_regular_q_transfer_v1 import sha

ROOT=Path(__file__).resolve().parent

def main():
    dest=ROOT/'shared_oblique_all12_inputs_v1'
    dest.mkdir()
    base=ROOT/'shared_oblique_collected_suffix_v1_mpi1'
    recover=ROOT/'shared_oblique_missing_fields_v1_mpi1'
    prior=ROOT/'shared_oblique_v_suffix_v1_mpi1'
    evidence=ROOT/'SHARED_OBLIQUE_MISSING_FIELDS_20261021.json'
    result=json.loads(evidence.read_text())
    assert result['independent_missing800_field_audit_pass'] and result['complete14_oblique_full_field_bindings']
    _,panels=wire(base/'unchanged_candidate_panels_rank0.bin')
    pm={(struct.unpack_from('<Q',b)[0],struct.unpack_from('<Q',b,32)[0]):b for b in panels}
    selection=[(recover/'recovered_full_field_points_rank0.bin',recover/'recovered_axis_records_rank0.bin'),
               (prior/'new_passed_full_field_points_rank0.bin',prior/'new_geometric_v_axis_proofs_rank0.bin'),
               (base/'new_passed_full_field_points_rank0.bin',base/'new_collected_axis_proofs_rank0.bin')]
    bindings={}; inputs=[evidence,base/'unchanged_candidate_panels_rank0.bin',
        base/'unchanged_frames_rank0.bin',base/'unchanged_actual_slabs_rank0.bin',
        ROOT/'SHARED_OBLIQUE_CLEARANCE_ADOPTION_20261021.json',
        ROOT/'SHARED_ATLAS_BUDGET_POLICY_20261021.json',
        ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin']
    for field_path,axis_path in selection:
        inputs.extend((field_path,axis_path))
        _,fields=wire(field_path);_,axes=wire(axis_path)
        am={struct.unpack_from('<3Q',b):i for i,b in enumerate(axes) if struct.unpack_from('<Q',b,24)[0]}
        for i,p in enumerate(fields):
            key=struct.unpack_from('<3Q',p);assert key not in bindings and key in am
            b=axes[am[key]]
            assert p[1424:1488]==b[256:320] and p[416:448]==b[224:256]
            bindings[key]=(field_path,i,axis_path,am[key])
    expected={(slot,index,axis) for slot,index in pm for axis in (0,1)}
    assert set(bindings)==expected and len(bindings)==8768
    per_slot=Counter(slot for slot,index in pm)
    with (dest/'complex_field_bindings.csv').open('w',newline='') as f:
        out=csv.writer(f)
        out.writerow(('slot_id','source','target','chart','panel','axis','full_field_packet','field_record','matching_axis_packet','axis_record'))
        for (slot,index,axis),link in sorted(bindings.items()):
            panel=pm[slot,index];source,target,chart=struct.unpack_from('<3Q',panel,8)
            out.writerow((slot,source,target,chart,index,axis,*map(str,link)))
    deps=dict(
        adopted14_subchart_geometry_density_clearance='complete',
        unchanged4384_typed_panel_binding='complete',
        all8768_complex_axis_full_field_radius_measure_inputs='complete',
        original_event_enclosures_retained='complete',
        shared_endpoint_error_attribution_to_all12='pending',
        real_center_moment_data_for_constant_removal='pending',
        explicit_frozen_CONTROL_thermal_context_binding='pending',
        general_evolved_heat_transport_and_restart='not_implemented',
        all12_analytic_remainders_and_native_forward_bounds='pending',
        complete_source_target_atlas_coverage='29_slots_pending',
        target_total_query_cost='unknown_for10_targets',
        native_R80_pairwise_sum_and_final_conversion='existing_implementation_not_executed_for_new_targets')
    contract=dict(scope='existing14 oblique subchart all12 INPUT dependency binding only; no new field or quadrature queries',
        complex_axis_field_binding_complete=True,oblique_subcharts=14,panels=4384,complex_axis_fields=8768,
        subchart_planned_positive_queries=4384*256,new_field_evaluations=0,new_physical_quadrature_queries=0,
        per_slot_panels=dict(sorted(per_slot.items())),dependencies=deps,
        remaining_incomplete_atlas_slots=29,all12_analytic_remainders_complete=False,
        full_quadrature_allowed=False,accepted_RK4_evolution=False,long_time_stability_demonstrated=False)
    (dest/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    (dest/'input.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in inputs))
    (dest/'source.sha256').write_text(sha(Path(__file__))+'  '+str(Path(__file__).resolve())+'\n')
    (dest/'output.sha256').write_text(''.join(sha(dest/n)+'  '+str(dest/n)+'\n' for n in ('contract.json','complex_field_bindings.csv')))
    print(json.dumps(contract,indent=2))

if __name__=='__main__':main()

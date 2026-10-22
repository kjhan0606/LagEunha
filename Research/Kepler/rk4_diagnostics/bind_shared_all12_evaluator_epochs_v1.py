"""Bind v6 packet IDs to their frozen evaluator manifests; no field calls."""
from pathlib import Path
import csv,json
from audit_shared_regular_q_transfer_v1 import rows,sha

ROOT=Path(__file__).resolve().parent

def main():
    registry=ROOT/'shared_atlas_adoption_v6';dest=ROOT/'shared_all12_evaluator_epochs_v1';dest.mkdir()
    packets=rows(registry/'typed_packet_registry.csv');inputs=[registry/'typed_packet_registry.csv',registry/'contract.json']
    state=ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin';state_sha=sha(state);assert state_sha=='74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005';inputs.append(state)
    helper_names=('material_shared_duffy.h','material_shared_oblique_point.h','material_shared_oblique_v_domain.h',
                  'material_shared_regular_prepared.h','material_shared_collected_radius.h',
                  'material_shared_weighted_angular.h','material_shared_origin_quadratic.h')
    result=[];manifests={}
    for p in packets:
        parent=Path(p['path']).parent
        if parent not in manifests:
            source=parent/'source.sha256';assert source.is_file()
            input_file=parent/'input.sha256'
            if not input_file.is_file():input_file=parent/'retained_input.sha256'
            assert input_file.is_file()
            binary=parent/'binary.sha256';assert binary.is_file()
            # Pin the native operation graph manifests, not just its output.
            value=dict(source_manifest=str(source),source_manifest_sha256=sha(source),
                input_manifest=str(input_file),input_manifest_sha256=sha(input_file),
                binary_manifest=str(binary),binary_manifest_sha256=sha(binary),
                original_state_sha256=state_sha,
                available_frozen_helper_sha256={name:sha(parent/name) for name in helper_names if (parent/name).is_file()})
            manifests[parent]=value;inputs.extend((source,input_file,binary))
        result.append(dict(packet_id=int(p['packet_id']),packet_sha256=p['sha256'],epoch=p['epoch'],**manifests[parent]))
    c=json.loads((registry/'contract.json').read_text())
    assert c['stored_new29_full_fields']==110332 and c['existing800_oblique_recovery_already_complete']
    contract=dict(scope='all v6 typed/field/proof packet IDs to original native source/input/binary closure manifests; helper presence is not a new usage or analytic proof',
        packet_evaluator_epochs=result,original_state_sha256=state_sha,
        preferred_new_missing_input_evaluator=str(ROOT/'shared_pending_clearance_v1_mpi1/source.sha256'),
        existing_saved_full_fields_to_reuse=110332,original800_oblique_fields_already_recovered=True,
        currently_unbound_complex_input_records=69910,
        list_requires_binding_existing_stronger_evidence_before_new_query=True,
        missing_source133_full_field_axis_records=7952,baseline_plan_axis_records_without_bound_complex_proof=61958,
        frozen_thermal_context_still_requires_explicit_binding=True,general_evolved_heat_restart_supported=False,
        new_field_evaluations=0,new_physical_quadrature_queries=0,
        all_target_analytic_certificates_complete=False,full_quadrature_allowed=False)
    (dest/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    (dest/'input.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in dict.fromkeys(inputs)))
    (dest/'source.py').write_bytes(Path(__file__).read_bytes())
    (dest/'source.sha256').write_text(sha(dest/'source.py')+'  '+str(dest/'source.py')+'\n')
    (dest/'output.sha256').write_text(sha(dest/'contract.json')+'  '+str(dest/'contract.json')+'\n')
    print(json.dumps(dict(bound_packet_epochs=len(result),bound_native_operation_graph_manifests=len(manifests),saved_full_fields=110332,unbound_complex_inputs=69910,new_queries=0)))

if __name__=='__main__':main()

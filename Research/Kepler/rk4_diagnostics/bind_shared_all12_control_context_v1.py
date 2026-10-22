"""Connect saved full fields to immutable CONTROL heat and source EOS contexts.

Only bounded metadata parsing. General evolved heat/restart is not supplied by
this frozen generating recipe and must never be inferred from a later current J.
"""
from pathlib import Path
from fractions import Fraction as F
import csv,json,struct
from audit_shared_regular_q_transfer_v1 import rows,sha,Q

ROOT=Path(__file__).resolve().parent
EXPECTED='74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005'

def main():
    dest=ROOT/'shared_all12_control_context_v1';dest.mkdir()
    registry=ROOT/'shared_atlas_adoption_v7'
    c=json.loads((registry/'contract.json').read_text())
    assert c['stored_active_full_fields']==180242 and c['unbound_or_missing_complex_axis_records']==0
    state=ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin'
    heat=ROOT/'shared_normalizers_v1_mpi4/immutable_heat_generating_state.bin'
    old_current=ROOT/'shared_normalizers_v1_mpi4/caller_unchanged_state.bin'
    for p in (state,heat,old_current):assert sha(p)==EXPECTED
    raw=state.read_bytes();h=struct.unpack_from('<8Q',raw)
    assert h[:3]==(0x4746535452433031,1,1) and h[4:7]==(384,1064,32)
    assert len(raw)==64+384*1064+h[3]*64
    call=ROOT/'shared_target_v3_mpi4'
    context=rows(call/'frozen_geometric_EOS_context.csv')
    original=rows(call/'source_target_contributions.csv')
    assert len(context)==24 and len(original)==384
    by_cell={int(r['cell']):r for r in context};assert len(by_cell)==24
    parameters=ROOT/'shared_density_witness_v1_mpi1/original_density_parameters.csv'
    param={int(r['parameter']):Q(r['value']) for r in rows(parameters)}
    assert param[1]>1 and param[2]>0
    for x in param.values():assert F(float(x))==x  # exact original binary64 parameter
    source_rows=[];index_to_cell={};seen=set()
    for r in original:
        s,cell=int(r['source']),int(r['cell']);assert s not in seen and 0<=s<384;seen.add(s)
        # STTriangle has the original native1064byte layout, final int cell.
        # Native load_saved in both audited new runs also verifies this index
        # and original packet physical-cell identity for every saved triangle.
        component=struct.unpack_from('<i',raw,64+1064*s+1060)[0]
        assert 0<=component<24 and cell in by_cell
        if component in index_to_cell:assert index_to_cell[component]==cell
        else:index_to_cell[component]=cell
        ctx=by_cell[cell];delta=float(ctx['passed_delta_double']);delta_dot=float(ctx['passed_delta_dot_double'])
        source_rows.append(dict(source=s,source_component=component,physical_source_cell=cell,
            generating_epoch=1,P0_hex=float(param[2]).hex(),gamma_hex=float(param[1]).hex(),
            delta_hex=delta.hex(),delta_dot_hex=delta_dot.hex(),
            delta_binary64_le=struct.pack('<d',delta).hex(),delta_dot_binary64_le=struct.pack('<d',delta_dot).hex(),
            immutable_generating_state_sha256=EXPECTED,current_CONTROL_state_sha256=EXPECTED))
    assert set(index_to_cell)==set(range(24)) and len(set(index_to_cell.values()))==24
    source_rows.sort(key=lambda r:r['source'])
    with (dest/'source_thermal_context_bindings.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(source_rows[0]));w.writeheader();w.writerows(source_rows)
    packets=rows(registry/'typed_packet_registry.csv');epochs={};bindings=[]
    helpers=('material_shared_regular_prepared.h','material_shared_collected_radius.h',
       'material_shared_weighted_angular.h','material_shared_origin_quadratic.h',
       'material_shared_oblique_point.h','material_shared_oblique_v_domain.h',
       'material_shared_duffy.h','material_shared_logistic_density.h')
    inputs=[registry/'contract.json',registry/'typed_packet_registry.csv',state,heat,old_current,
            call/'source_target_contributions.csv',call/'frozen_geometric_EOS_context.csv',parameters]
    for p in packets:
        parent=Path(p['path']).parent
        if parent not in epochs:
            source=parent/'source.sha256';native_input=parent/'input.sha256'
            if not native_input.is_file():native_input=parent/'retained_input.sha256'
            binary=parent/'binary.sha256'
            assert source.is_file() and native_input.is_file() and binary.is_file()
            epochs[parent]=dict(source_manifest=str(source),source_manifest_sha256=sha(source),
                input_manifest=str(native_input),input_manifest_sha256=sha(native_input),
                binary_manifest=str(binary),binary_manifest_sha256=sha(binary),
                available_frozen_helper_sha256={name:sha(parent/name) for name in helpers if (parent/name).is_file()})
            inputs.extend((source,native_input,binary))
        bindings.append(dict(packet_id=int(p['packet_id']),packet_sha256=p['sha256'],epoch=p['epoch'],**epochs[parent]))
    contract=dict(scope='all v7 packets to frozen native epochs and all384source EOS/CONTROL heat contexts; metadata/identity binding only, no new field or moment queries',
       active_full_fields_available=180242,packet_evaluator_epochs=bindings,
       frozen_CONTROL_thermal_context_bound=True,source_context_rows=384,source_cells=24,
       current_state_sha256=EXPECTED,immutable_heat_generating_state_sha256=EXPECTED,
       physical_reference_heat='P0/(gamma-1) * J_generating^(1-gamma)',
       use_current_J_for_generating_heat_allowed_only_for_this_byte_identical_CONTROL_state=True,
       source_EOS_delta_and_delta_dot='original retained binary64 values, indexed by source component, not recipient target',
       helper_presence_is_not_proof_of_usage=True,existing_native_load_saved_validates_source_cell_identity=True,
       required_reference_density_evaluator=c['required_reference_density_evaluator'],
       required_reference_density_evaluator_sha256=c['required_reference_density_evaluator_sha256'],
       general_evolved_heat_restart_supported=False,heat_reset_permitted=False,
       new_complex_field_queries=0,new_moment_queries=0,new_physical_quadrature_queries=0,
       all_target_analytic_certificates_complete=False,full_quadrature_allowed=False,
       accepted_RK4_evolution=False,long_Kepler_validated=False)
    (dest/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    (dest/'input.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in dict.fromkeys(inputs)))
    (dest/'source.py').write_bytes(Path(__file__).read_bytes())
    (dest/'source.sha256').write_text(sha(dest/'source.py')+'  '+str(dest/'source.py')+'\n')
    (dest/'output.sha256').write_text(''.join(sha(dest/n)+'  '+str(dest/n)+'\n' for n in ('contract.json','source_thermal_context_bindings.csv')))
    print(json.dumps(dict(bound_source_contexts=384,bound_packet_epochs=len(bindings),native_epochs=len(epochs),
        available_fields=180242,frozen_CONTROL_thermal_context_bound=True,
        general_evolved_heat_restart_supported=False,new_field_or_moment_queries=0)))

if __name__=='__main__':main()

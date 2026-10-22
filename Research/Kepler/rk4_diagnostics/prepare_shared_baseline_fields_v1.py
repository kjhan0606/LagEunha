"""Bounded inventory of the v6 missing-field epochs, then an exact recovery list.

No directory traversal, geometry evaluation or physical quadrature. Existing
source133 clearance certificates are retained; these jobs recover absent data.
"""
from collections import Counter
from pathlib import Path
import csv, json, struct
from audit_shared_regular_q_transfer_v1 import rows, sha

ROOT = Path(__file__).resolve().parent

def main():
    dest = ROOT/'shared_baseline_field_inputs_v1'
    dest.mkdir()
    registry = ROOT/'shared_atlas_adoption_v6'
    inputs = [registry/n for n in ('contract.json', 'effective_slot_plan_references.csv', 'typed_packet_registry.csv')]
    packets = {p['path']:p for p in rows(registry/'typed_packet_registry.csv')}
    selected = [r for r in rows(registry/'effective_slot_plan_references.csv')
                if r['full_field_evidence_status'] in ('unbound','missing_snapshot')]
    paths = sorted({r['typed_packet'] for r in selected})
    inventory = []
    for parent in sorted({str(Path(p).parent) for p in paths}):
        folder = Path(parent)
        manifest = folder/'output.sha256'
        inputs.append(manifest)
        outputs = []
        for line in manifest.read_text().splitlines():
            digest, name = line.split('  ', 1)
            outputs.append(Path(name).name)
        # The six enumerated native producer epochs only wrote plans, state,
        # scalar clearance summaries and the narrow source133 axis proofs.
        # In particular none wrote an sd::Point packet. Unknown output names
        # require inspection before they could be silently classified missing.
        expected = {
          'shared_all_atlas_v2_mpi4': {'contract.json','all_target_cost.csv','caller_unchanged_state.bin'} |
              {f'{stem}_rank{rank}.{ext}' for rank in range(4) for stem,ext in
               (('actual_atlas_cost','csv'),('actual_chart_status','csv'),('typed_new_atlas','bin'))},
          'shared_origin_atlas_v2_mpi1': {'contract.json','origin_atlas_cost_and_clearance.csv','graded_source133_atlas_rank0.bin'},
          'shared_origin_atlas_v3_mpi1': {'contract.json','origin_atlas_cost_and_clearance.csv','graded_source133_atlas_rank0.bin'},
          'shared_atlas_repairs_v1_mpi1': {'contract.json','cost_repairs.csv','oblique_direction_proofs.csv','vertex_seed_inverse.csv','cost_repaired_panels_rank0.bin','vertex_seed_inverse_rank0.bin'},
          'shared_ray_suffix_v1_mpi1': {'contract.json','first_failure_arguments.csv','suffix_clearance.csv','suffix_axis_proofs_rank0.bin','unchanged_original_candidate_plans_rank0.bin'},
          'shared_weighted_ray_v1_mpi1': {'contract.json','weighted_axis_clearance.csv','new_weighted_axis_proofs_rank0.bin','unchanged_original_plan_rank0.bin','unchanged_previous_passing_axes_rank0.bin'},
        }
        assert set(outputs) == expected[folder.name], (folder, outputs)
        inventory.append(dict(epoch=parent, output_manifest_sha256=sha(manifest),
                              outputs=outputs, full_field_packet_written=False,
                              existing_scalar_or_axis_clearance_preserved=True))
    cache = {}
    for path in paths:
        data=Path(path).read_bytes(); h=struct.unpack_from('<8Q',data)
        assert h[3]==704 and len(data)==64+704*h[4] and sha(Path(path))==packets[path]['sha256']
        cache[path]=data; inputs.append(Path(path))
    jobs=[]; counts=Counter()
    for r in selected:
        first,n=int(r['typed_record_first']),int(r['typed_record_count'])
        slot,source,target=int(r['slot_id']),int(r['source']),int(r['target'])
        assert n==int(r['active_panels']) and n>0
        for j in range(n):
            actual=struct.unpack_from('<4Q',cache[r['typed_packet']],64+704*(first+j))
            assert actual==(slot,source,j,target), actual
        kind='recover_historical_pass_data' if r['full_field_evidence_status']=='missing_snapshot' else 'new_complex_evidence_and_full_field'
        counts[kind]+=2*n
        jobs.append(dict(slot_id=slot,source=source,target=target,typed_packet=r['typed_packet'],
                         record_first=first,record_count=n,work_kind=kind))
    assert sum(2*j['record_count'] for j in jobs)==69910
    assert counts=={'new_complex_evidence_and_full_field':61958,'recover_historical_pass_data':7952}
    with (dest/'slot_jobs.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(jobs[0]));w.writeheader();w.writerows(jobs)
    contract=dict(scope='enumerated six v6 baseline native producer epochs; no saved full-field packets in their output manifests; exact missing-field recovery inputs, not replay of the existing110332 fields',
       producer_inventory=inventory,slot_jobs=len(jobs),regular_panels=34955,
       maximum_complex_axis_queries=69910,axis_work_kind=dict(counts),
       maximum_full_field_packet_bytes=64+69910*1680,
       maximum_axis_proof_packet_bytes=64+69910*624,
       maximum_boundary_packet_bytes=64+69910*1552,
       maximum_output_bytes_conservative=400_000_000,
       fixed_original_plans=True,selected_branch_gradient_cones_only=True,
       new_arrangements=0,new_event_root_jobs=0,new_inverse_solves=0,
       original800_oblique_fields_already_available=True,new_physical_quadrature_queries=0,
       frozen_CONTROL_state_only=True,general_evolved_heat_restart_supported=False,
       full_quadrature_allowed=False,long_Kepler_validated=False)
    (dest/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    (dest/'input.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in dict.fromkeys(inputs)))
    (dest/'source.py').write_bytes(Path(__file__).read_bytes())
    (dest/'source.sha256').write_text(sha(dest/'source.py')+'  '+str(dest/'source.py')+'\n')
    (dest/'output.sha256').write_text(''.join(sha(dest/n)+'  '+str(dest/n)+'\n' for n in ('slot_jobs.csv','contract.json')))
    print(json.dumps({k:v for k,v in contract.items() if k!='producer_inventory'},indent=2))

if __name__=='__main__':main()

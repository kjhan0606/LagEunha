"""Adopt completed readonly target integrals, preserving the remaining transfer."""
from pathlib import Path
import csv,json,struct
from audit_shared_regular_q_transfer_v1 import rows,sha

ROOT=Path(__file__).resolve().parent

def main():
    dest=ROOT/'shared_native_transfer_adoption_v1';dest.mkdir()
    run=ROOT/'shared_whole_native_v1_mpi1'
    audit=ROOT/'shared_whole_native_audit_v1/contract.json'
    a=json.loads(audit.read_text());assert a['independent_whole_native_saved_accounting_pass'] and a['readonly_target_candidates_pass']==4
    registry=ROOT/'shared_atlas_adoption_v7'
    slots=rows(registry/'effective_slot_plan_references.csv')
    active={int(r['slot_id']):r for r in slots if int(r['active_panels'])}
    data=(run/'native_source_packets_rank0.bin').read_bytes();h=struct.unpack_from('<8Q',data)
    assert h[3:5]==(832,64)
    sources=[];selected=set(a['selected_targets']);source_ids=set()
    for j in range(64):
        slot,source,target,planned,count,complete,oracle=struct.unpack_from('<7Q',data,64+832*j)
        assert complete==1 and planned==count==22528 and target in selected
        assert int(active[slot]['source'])==source and int(active[slot]['target'])==target
        # These sources have one complete whole recipient and no other active
        # target overlap. This is a64source fact, not full384source closure.
        assert [int(r['slot_id']) for r in slots if int(r['source'])==source and int(r['active_panels'])]==[slot]
        source_ids.add(source)
        sources.append(dict(slot_id=slot,source=source,target=target,native_packet=str(run/'native_source_packets_rank0.bin'),native_record=j,
           positive_packet=str(run/f'positive_slot{slot}.bin'),points=count,
           outgoing_incoming_overlap_record='same_slot_same_native_packet_record',
           whole_source_coverage=True))
    assert len(source_ids)==64
    completed=[]
    for target in sorted(selected):
        path=run/f'readonly_target_packet_{target}.bin'
        assert struct.unpack_from('<2Q',path.read_bytes())==(1,target)
        completed.append(dict(target=target,path=str(path),sha256=sha(path),status='readonly_original_total_budget_and_positivity_pass',queries=360448))
    original=ROOT/'shared_missing_source_v1_mpi1/caller_target_packet.bin'
    assert sha(original)=='946b5639159c57f34d9104a7aac537f6dcfa6a0e232d40057ac5aa465df532f3'
    completed.append(dict(target=68526,path=str(original),sha256=sha(original),status='retained_original_whole_target_publication_not_split_into_sources',queries=268032))
    remaining=sorted({int(r['target']) for r in slots}-{r['target'] for r in completed})
    remaining_slots=[r for r in slots if int(r['target']) in remaining and int(r['active_panels'])]
    queries=sum(int(r['active_panels'])*256 for r in remaining_slots)
    assert len(remaining)==19 and queries==21629184
    for name,records in (('completed_target_references.csv',completed),('completed_source_overlap_references.csv',sources)):
        with (dest/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');w.writeheader();w.writerows(records)
    report=dict(scope='four audited readonly whole-target native integrals plus original68526publication; actual overlap values/errors bound for64single-recipient sources; no full transfer or physical state commit',
       actual_targets_with_complete_original_error_and_positivity_evidence=5,
       new_readonly_targets=sorted(selected),original68526publication_preserved=True,
       completed_native_source_overlap_records=64,source_single_recipient_whole_coverage_verified=64,
       newly_acquired_actual_native_queries=1441792,retained_original_queries=268032,
       remaining_target_ids=remaining,remaining_nonempty_slots=len(remaining_slots),
       remaining_planned_actual_queries=queries,
       maximum_new_whole_target_total_error_fraction=a['worst_target_budgets'][0],
       completed_actual_results_must_not_be_replayed=True,
       source_outgoing_target_incoming_reference_same_overlap_packet=True,
       original_whole4_group_not_divided=True,
       full384source_numerical_closure_complete=False,full24target_transfer_complete=False,
       original_caller_state_committed=False,general_evolved_heat_restart_supported=False,
       accepted_RK4_evolution=False,long_Kepler_validated=False)
    (dest/'contract.json').write_text(json.dumps(report,indent=2)+'\n')
    inputs=[audit,registry/'contract.json',registry/'effective_slot_plan_references.csv',run/'native_source_packets_rank0.bin',run/'output.sha256',original,
            *(Path(r['path']) for r in completed)]
    (dest/'input.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in dict.fromkeys(inputs)))
    (dest/'source.py').write_bytes(Path(__file__).read_bytes())
    (dest/'source.sha256').write_text(sha(dest/'source.py')+'  '+str(dest/'source.py')+'\n')
    (dest/'output.sha256').write_text(''.join(sha(dest/n)+'  '+str(dest/n)+'\n' for n in ('contract.json','completed_target_references.csv','completed_source_overlap_references.csv')))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()

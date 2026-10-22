"""Select complete whole-only targets for bounded native error acquisition.

The rest of the24-target transfer remains required. No physical or field query.
"""
from pathlib import Path
from collections import defaultdict
import csv,json,struct
from audit_shared_regular_q_transfer_v1 import rows,sha

ROOT=Path(__file__).resolve().parent

def main():
    dest=ROOT/'shared_whole_native_inputs_v1';dest.mkdir()
    registry=ROOT/'shared_atlas_adoption_v7'
    evidence=ROOT/'shared_saved_all12_evidence_v1/contract.json'
    e=json.loads(evidence.read_text());assert e['all_new23_target_all12_analytic_subtotals_within_original_tau']
    domain=ROOT/'SHARED_REGISTRY_DOMAINS_BRIDGED_20261021.json'
    pairs=ROOT/'SHARED_SOURCE_PAIRS_20261021.json'
    d=json.loads(domain.read_text());g=json.loads(pairs.read_text())
    assert d['effective_all9216_domain_arrangements_complete'] and d['saved_event_source_partition_and_edge_cancellation_complete']
    assert g['global_source_image_interiors_disjoint_complete']
    packets={int(r['packet_id']):r for r in rows(registry/'typed_packet_registry.csv')}
    raw={};by_target=defaultdict(list);mixed=set()
    inputs=[evidence,domain,pairs,registry/'contract.json',registry/'panel_evidence_registry.csv',registry/'typed_packet_registry.csv']
    for r in rows(registry/'panel_evidence_registry.csv'):
        target=int(r['target'])
        if target==68526:continue
        by_target[target].append(r)
        if r['panel_kind']!='regular':mixed.add(target);continue
        pid=int(r['typed_packet_id'])
        if pid not in raw:
            p=Path(packets[pid]['path']);data=p.read_bytes();h=struct.unpack_from('<8Q',data)
            assert h[3]==704 and len(data)==64+704*h[4] and sha(p)==packets[pid]['sha256']
            raw[pid]=data;inputs.append(p)
        first=64+704*int(r['typed_record'])
        if raw[pid][first+192]!=1:mixed.add(target)
    selected=sorted(set(by_target)-mixed);assert selected
    jobs=[];work=[]
    slot_rows=rows(registry/'effective_slot_plan_references.csv')
    for target in selected:
        n=len(by_target[target]);cost=n*256;assert cost<=13074432
        work.append(dict(target=target,panels=n,planned_positive_queries=cost,maximum_positive_packet_bytes=1008*cost))
        for r in slot_rows:
            if int(r['target'])==target and int(r['active_panels']):
                assert r['full_field_evidence_status']=='stored_pair'
                jobs.append({k:r[k] for k in ('slot_id','source','target','typed_packet','typed_record_first','typed_record_count')})
    with (dest/'slot_jobs.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(jobs[0]),lineterminator='\n');w.writeheader();w.writerows(jobs)
    report=dict(scope='complete whole-only target read-only native numerical-error acquisition; all other targets and full conservative transfer remain required',
       selected_targets=work,selected_slots=len(jobs),selected_panels=sum(r['panels'] for r in work),
       maximum_new_native_queries=sum(r['planned_positive_queries'] for r in work),
       maximum_positive_packet_bytes=sum(r['maximum_positive_packet_bytes'] for r in work),
       input_field_replays=0,old68526publication_replayed=False,
       source_partition_and_pair_geometry_evidence_reused=True,
       actual_native_error_unknown_until_acquired=True,
       new_physical_queries_in_this_preparation=0,full24target_quadrature_allowed=False,
       caller_state_commit_allowed=False,accepted_RK4_evolution=False,long_Kepler_validated=False)
    (dest/'contract.json').write_text(json.dumps(report,indent=2)+'\n')
    (dest/'input.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in dict.fromkeys(inputs)))
    (dest/'source.py').write_bytes(Path(__file__).read_bytes())
    (dest/'source.sha256').write_text(sha(dest/'source.py')+'  '+str(dest/'source.py')+'\n')
    (dest/'output.sha256').write_text(''.join(sha(dest/n)+'  '+str(dest/n)+'\n' for n in ('contract.json','slot_jobs.csv')))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()

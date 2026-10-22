"""Assemble all existing overlap packets; no fields, refit or state mutation.

Source partition integrals are obtained from certified complete overlaps.
This is not an independent whole-source integration or an evolved state.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import csv,json,struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_regular_q_transfer_v1 import rows,Q,extended,sha

ROOT=Path(__file__).resolve().parent
NAMES=('mass','px','py','angular','kinetic','physical_U','volume','qx','qy','polar','physical_Udot','potential')

def write_csv(path,records):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');w.writeheader();w.writerows(records)

def main():
    dest=ROOT/'shared_transfer_incidence_v1';dest.mkdir()
    whole=ROOT/'shared_whole_native_v1_mpi1';mixed=ROOT/'shared_mixed_native_v2_mpi1'
    audit=ROOT/'shared_mixed_native_summary_audit_v1/contract.json'
    c=json.loads(audit.read_text());assert c['saved_mixed_runtime_accounting_pass'] and c['all_remaining19_total_budgets_pass'] and not c['failed_oracle_comparisons'] and c['remaining_unacquired_queries']==0
    oldaudit=ROOT/'SHARED_MISSING_SOURCE_LEDGER_20261021.json';old=json.loads(oldaudit.read_text());assert old['target_published'] and not old['ORIGINAL_tau_failures']
    recovered=ROOT/'shared_original_whole4_source_accounting_v1'
    recovery=json.loads((recovered/'contract.json').read_text());assert recovery['all12_original_group_reconciliations_pass']
    slots=rows(ROOT/'shared_atlas_adoption_v7/effective_slot_plan_references.csv')
    active={int(r['slot_id']):r for r in slots if int(r['active_panels'])};assert len(slots)==9216 and len(active)==559
    targets=sorted({int(r['target']) for r in slots});assert len(targets)==24
    overlap={};inputs=[audit,oldaudit,recovered/'contract.json',recovered/'original_whole4_source_accounting.csv',ROOT/'shared_atlas_adoption_v7/effective_slot_plan_references.csv',ROOT/'SHARED_REGISTRY_DOMAINS_BRIDGED_20261021.json',ROOT/'SHARED_SOURCE_PAIRS_20261021.json']
    def add(slot,source,target,k,n,e,remainder,points,provenance):
        assert slot==24*source+targets.index(target) and slot in active and points==256*int(active[slot]['active_panels'])
        assert e>=0 and remainder>=0
        if slot not in overlap:overlap[slot]=dict(source=source,target=target,points=points,provenance=provenance,values={})
        r=overlap[slot];assert (r['source'],r['target'],r['points'])==(source,target,points) and k not in r['values']
        r['values'][k]=(n,e,remainder)
    for run in (whole,mixed):
        path=run/'native_source_packets_rank0.bin';inputs.extend([path,run/'output.sha256'])
        _,packets=wire(path)
        for record,p in enumerate(packets):
            slot,source,target,planned,points,complete,oracles=struct.unpack_from('<7Q',p)
            assert len(p)==832 and complete and points==planned
            for k in range(12):add(slot,source,target,k,quad(p,64+32*k),quad(p,80+32*k),quad(p,448+16*k),points,f'{path}#{record}')
    oldrun=ROOT/'shared_missing_source_v1_mpi1';path=oldrun/'actual_clipped_source_integrals.csv';inputs.append(path)
    clipped=rows(path);assert len(clipped)==26*12
    for record,r in enumerate(clipped):
        source=int(r['source']);assert r['all12_complete']=='1'
        add(24*source+targets.index(68526),source,68526,NAMES.index(r['moment']),Q(r['value']),Q(r['native_geometry_pairwise_sum_error']),Q(r['analytic_remainder']),int(r['positive_points']),f'{path}#row{record}')
    path=recovered/'original_whole4_source_accounting.csv'
    for record,r in enumerate(rows(path)):
        add(int(r['slot_id']),int(r['source']),68526,NAMES.index(r['moment']),F(r['exact_saved_native_center']),F(r['exact_sum_point_radii']),F(r['exact_sum_original_analytic_remainders']),int(r['saved_positive_points']),f'{path}#row{record}')
    assert set(overlap)==set(active) and all(set(r['values'])==set(range(12)) for r in overlap.values())
    assert sum(r['points'] for r in overlap.values())==23339008
    outgoing=defaultdict(lambda:[F(0) for _ in NAMES]);incoming=defaultdict(lambda:[F(0) for _ in NAMES]);source_error=defaultdict(lambda:[F(0) for _ in NAMES]);target_error=defaultdict(lambda:[F(0) for _ in NAMES]);overlap_rows=[]
    for slot,r in sorted(overlap.items()):
        for k,name in enumerate(NAMES):
            n,e,remainder=r['values'][k];s,t=r['source'],r['target'];outgoing[s][k]+=n;incoming[t][k]+=n;source_error[s][k]+=e+remainder;target_error[t][k]+=e+remainder
            overlap_rows.append(dict(slot_id=slot,source=s,target=t,moment=name,points=r['points'],exact_center=str(n),exact_native_radius=str(e),exact_analytic_remainder=str(remainder),outgoing_incoming_same_record=True,provenance=r['provenance']))
    assert set(outgoing)==set(range(384)) and set(incoming)==set(targets)
    target_rows=[];global_report=[];published_centers=defaultdict(lambda:[F(0) for _ in NAMES]);published_errors=defaultdict(lambda:[F(0) for _ in NAMES])
    for t in targets:
        path=oldrun/'caller_target_packet.bin' if t==68526 else (whole if t in (17894,20502,39878,44945) else mixed)/f'readonly_target_packet_{t}.bin'
        inputs.append(path);data=path.read_bytes();assert len(data)==400 and struct.unpack_from('<2Q',data)==(1,t)
        if t==68526:assert sha(path)=='946b5639159c57f34d9104a7aac537f6dcfa6a0e232d40057ac5aa465df532f3'
        for k in range(12):
            n=extended(data[16+16*k:32+16*k]);e=extended(data[208+16*k:224+16*k]);assert abs(n-incoming[t][k])+target_error[t][k]<=e
            if k in (0,5):assert n>e
            published_centers[t][k]=n;published_errors[t][k]=e
        target_rows.append(dict(target=t,original_packet=str(path),sha256=sha(path),state='complete_original_budget_readonly_packet',overlap_count=sum(r['target']==t for r in overlap.values())))
    source_rows=[]
    for s in range(384):
        for k,name in enumerate(NAMES):
            if k in (0,5):assert outgoing[s][k]>source_error[s][k]
            source_rows.append(dict(source=s,moment=name,exact_partition_integral_center=str(outgoing[s][k]),exact_native_plus_analytic_radius=str(source_error[s][k]),overlap_count=sum(r['source']==s for r in overlap.values()),independent_whole_source_comparison=False))
    for k,name in enumerate(NAMES):
        before=sum(outgoing[s][k] for s in range(384));after=sum(incoming[t][k] for t in targets);assert before==after
        rounded=sum(published_centers[t][k] for t in targets);error=sum(published_errors[t][k] for t in targets)
        assert abs(rounded-before)+sum(target_error[t][k] for t in targets)<=error
        global_report.append(dict(moment=name,exact_overlap_outgoing_incoming_difference='0',exact_target_serialization_difference=str(rounded-before),total_published_error=str(error)))
    write_csv(dest/'all559_overlap_moments.csv',overlap_rows);write_csv(dest/'all384_source_partition_integrals.csv',source_rows);write_csv(dest/'all24_completed_target_references.csv',target_rows)
    result=dict(scope='all559original overlaps,384source partition integrals and24original-budget target packets linked; exact shared-overlap incidence and final target serialization accounting; no independent whole-source before/after comparison, no state evolution or refit',
        all24_original_target_budget_packets_complete=True,nonempty_overlap_packets=559,source_partition_integrals=384,
        actual_saved_queries_reused=23339008,new_field_or_physical_queries=0,
        original68526packet_preserved=True,outgoing_incoming_same_overlap_packet=True,
        exact_all12_shared_overlap_incidence_pass=True,all384_source_partition_mass_and_U_positive=True,
        all24_target_serialization_enclosures_pass=True,global_accounting=global_report,
        full384_independent_whole_source_comparison_complete=False,full24target_transfer_state_committed=False,
        retained_field_fragments_implemented=False,general_evolved_heat_restart_supported=False,
        original_caller_state_committed=False,accepted_RK4_evolution=False,long_Kepler_validated=False)
    (dest/'contract.json').write_text(json.dumps(result,indent=2)+'\n');(dest/'source.py').write_bytes(Path(__file__).read_bytes())
    (dest/'source.sha256').write_text(sha(dest/'source.py')+'  '+str(dest/'source.py')+'\n')
    (dest/'input.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in dict.fromkeys(inputs)))
    (dest/'output.sha256').write_text(''.join(sha(dest/n)+'  '+str(dest/n)+'\n' for n in ('contract.json','all559_overlap_moments.csv','all384_source_partition_integrals.csv','all24_completed_target_references.csv')))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()

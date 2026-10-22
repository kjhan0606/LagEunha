"""Compose immutable original/replacement plan references; no geometry replay."""
from collections import Counter,defaultdict
from pathlib import Path
import csv
import json
import struct
import sys
from audit_shared_regular_q_transfer_v1 import rows,sha

ROOT=Path(__file__).resolve().parent
def write_csv(path,records):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)

def main():
    dest=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_atlas_adoption_v1')
    assert dest.parent==ROOT and not dest.exists()
    auditpaths=[ROOT/name for name in ('SHARED_ALL_ATLAS_COST_20261021.json','SHARED_ATLAS_WITNESS_20261021.json','SHARED_ATLAS_REPAIRS_20261021.json','SHARED_NORMALIZERS_20261021.json')]
    audits=[json.loads(p.read_text()) for p in auditpaths]
    assert audits[0]['packet_and_cost_accounting_pass'] and audits[1]['independent_saved_witness_accounting_pass'] and audits[2]['limited_saved_repair_audit_pass']
    review=ROOT/'ASTRA_SHARED_FAILED_ATLAS_END_REVIEW_20261021.md';assert review.is_file()
    old=ROOT/'shared_all_atlas_v2_mpi4';new=ROOT/'shared_atlas_repairs_v1_mpi1'
    for run in (old,new):
        for manifest in ('source.sha256','output.sha256'):
            for line in (run/manifest).read_text().splitlines():
                digest,path=line.split('  ',1);assert sha(Path(path))==digest
    original={int(r['slot_id']):r for k in range(4) for r in rows(old/f'actual_atlas_cost_rank{k}.csv')};assert len(original)==9216
    repairs={int(r['slot_id']):r for r in rows(new/'cost_repairs.csv')};assert set(repairs)=={3248,4499}
    ranges={};artifacts=[]
    for path in [*(old/f'typed_new_atlas_rank{k}.bin' for k in range(4)),new/'cost_repaired_panels_rank0.bin']:
        raw=path.read_bytes();h=struct.unpack_from('<8Q',raw);assert h[3]==704 and len(raw)==64+h[4]*704
        artifacts.append({'path':str(path),'SHA256':sha(path),'bytes':len(raw),'record_count':h[4]})
        for j in range(h[4]):
            slot,source,index,target=struct.unpack_from('<4Q',raw,64+j*704)
            key=str(path),slot
            if key not in ranges:assert index==0;ranges[key]=[j,0]
            first,count=ranges[key];assert index==count and j==first+count
            assert source==slot//24 and target==int(original[slot]['target']);ranges[key][1]+=1
    records=[];cost=defaultdict(Counter)
    for slot,r in sorted(original.items()):
        replacement=slot in repairs;reused=r['reused']=='1';complete=replacement or r['atlas_complete']=='1'
        count=int(repairs[slot]['panels'] if replacement else r['stored_panels'])
        if replacement:assert r['atlas_complete']=='0' and int(r['stored_panels'])==3192 and repairs[slot]['slot_cost_pass']=='1'
        kind='accepted_replacement_cost_plan' if replacement else 'retained_original_target_plan' if reused else 'retained_original_other_target_plan'
        path=new/'cost_repaired_panels_rank0.bin' if replacement else old/f'typed_new_atlas_rank{slot%4}.bin'
        first=-1
        if count and not reused:
            first,savedcount=ranges[str(path),slot];assert savedcount==count
        target=int(r['target']);v=cost[target];v['incomplete']+=not complete;v['retained_panels' if reused else 'new_panels']+=count
        records.append({'slot_id':slot,'source':int(r['source']),'target':target,'atlas_complete':int(complete),'active_panels':count,'complete_planned_queries':count*256 if complete else 'UNKNOWN','provenance':kind,'typed_packet':str(path) if count and not reused else 'ORIGINAL_TARGET_ATLAS_REFERENCE' if reused and count else 'NONE','typed_record_first':first,'typed_record_count':count if not reused else 0,'original_atlas_complete':r['atlas_complete'],'original_stored_panels':r['stored_panels'],'original_reason':r['reason'],'analytic_ready':r['analytic_ready']})
    targets=[]
    for target,v in sorted(cost.items()):
        planned=(v['new_panels']+v['retained_panels'])*256;complete=v['incomplete']==0
        targets.append({'target':target,'incomplete_atlas_slots':v['incomplete'],'active_new_panels':v['new_panels'],'retained_original_target_panels':v['retained_panels'],'complete_planned_queries':planned if complete else 'UNKNOWN','query_cap':817152,'complete_cost_pass':int(complete and planned<=817152),'actual_new_query_packet_bytes':v['new_panels']*256*1008 if complete else 'UNKNOWN','all_target_analytic_certificates_complete':0})
    assert sum(r['incomplete_atlas_slots'] for r in targets)==33
    assert sum(r['active_new_panels'] for r in targets)==30979
    assert sum(r['retained_original_target_panels'] for r in targets)==1047
    assert all(r['complete_cost_pass']==1 for r in targets if r['incomplete_atlas_slots']==0)
    dest.mkdir();write_csv(dest/'effective_slot_plan_references.csv',records);write_csv(dest/'effective_target_cost.csv',targets)
    paths=[Path(__file__).resolve(),*auditpaths,review,*(Path(a['path']) for a in artifacts),*(old/f'actual_atlas_cost_rank{k}.csv' for k in range(4)),new/'cost_repairs.csv']
    (dest/'input.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in paths))
    result={'scope':'immutable plan-reference and exact cost composition after accepted2 replacements, original partial prefixes remain archived','effective_atlas_slots':9216,'incomplete_atlas_slots':33,'effective_active_new_panels':30979,'retained_original_target_panels':1047,'archived_replaced_partial_panels':6384,'new_replacement_panels':256,'normalizers_other23_complete_from_separate_audit':True,'original68526_tau_and_packet_unchanged':True,'complete_target_costs_within_original_budget':True,'incomplete_target_costs_unknown':True,'all_target_actual_atlas_complete':False,'all_target_analytic_certificates_complete':False,'full_quadrature_allowed':False,'new_physical_quadrature_queries':0,'state_or_heat_changes':0,'artifacts':artifacts,'target_cost':targets,'review_path':str(review),'review_SHA256':sha(review)}
    (dest/'contract.json').write_text(json.dumps(result,indent=2)+'\n')
    (dest/'output.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in (dest/'contract.json',dest/'effective_slot_plan_references.csv',dest/'effective_target_cost.csv')))
    print(json.dumps({k:result[k] for k in ('scope','incomplete_atlas_slots','effective_active_new_panels','retained_original_target_panels','archived_replaced_partial_panels','new_replacement_panels','complete_target_costs_within_original_budget','incomplete_target_costs_unknown','full_quadrature_allowed')},indent=2))

if __name__=='__main__':main()

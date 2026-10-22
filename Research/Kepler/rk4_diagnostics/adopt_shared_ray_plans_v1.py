"""Adopt only independently covered and clearance-complete same-ray plans."""
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
    dest=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_atlas_adoption_v4');assert dest.parent==ROOT and not dest.exists()
    previous=ROOT/'shared_atlas_adoption_v3';coverage=ROOT/'shared_ray_coverage_v1';run=ROOT/'shared_ray_suffix_v1_mpi1'
    for directory in (previous,coverage,run):
        for name in ('input.sha256','output.sha256'):
            for line in (directory/name).read_text().splitlines():
                digest,path=line.split('  ',1);assert sha(Path(path))==digest
    adopted=json.loads((coverage/'contract.json').read_text());assert adopted['domain_coverage_and_native_geometry_density_clearance_complete']
    ready={r['slot_id']:r for r in adopted['adoptable_slot_plans']};assert set(ready)=={3200}
    policy=ROOT/'SHARED_ATLAS_BUDGET_POLICY_20261021.json';p=json.loads(policy.read_text());limit=p['new_per_target_query_limit'];assert limit==13074432
    path=run/'unchanged_original_candidate_plans_rank0.bin';raw=path.read_bytes();h=struct.unpack_from('<8Q',raw);assert h[3]==704
    ranges={}
    for j in range(h[4]):
        slot,source,index,target=struct.unpack_from('<4Q',raw,64+j*704)
        if slot not in ready:continue
        if slot not in ranges:assert index==0;ranges[slot]=[j,0]
        first,count=ranges[slot];assert index==count and j==first+count;ranges[slot][1]+=1
    ledger=rows(previous/'effective_slot_plan_references.csv');cost=defaultdict(Counter);pending=[]
    for r in ledger:
        slot=int(r['slot_id']);r['required_reference_radius_evaluator']='original_reference_radius'
        if slot in ready:
            proof=ready[slot];first,count=ranges[slot];assert count==proof['panels'] and r['atlas_complete']=='0' and r['active_panels']=='0'
            r.update(atlas_complete='1',active_panels=str(count),complete_planned_queries=str(count*256),provenance='adopted_actual_plan_same_function_collected_ray_clearance',typed_packet=str(path),typed_record_first=str(first),typed_record_count=str(count),required_reference_radius_evaluator='collected_reference_radius_same_exact_function')
        v=cost[int(r['target'])];v['pending']+=r['atlas_complete']=='0';v['retained' if r['provenance']=='retained_original_target_plan' else 'new']+=int(r['active_panels'])
        if r['atlas_complete']=='0':pending.append({'slot_id':slot,'source':int(r['source']),'target':int(r['target']),'original_failure':r['original_reason']})
    targets=[]
    for target,v in sorted(cost.items()):
        known=(v['new']+v['retained'])*256;targets.append({'target':target,'incomplete_atlas_slots':v['pending'],'active_new_panels':v['new'],'retained_original_target_panels':v['retained'],'known_positive_query_lower_bound':known,'complete_planned_queries':'UNKNOWN' if v['pending'] else known,'authorized_query_cap':limit,'complete_cost_pass':int(not v['pending'] and known<=limit),'analytic_certificates_complete':0})
    assert len(ledger)==9216 and len(pending)==30 and sum(r['active_new_panels'] for r in targets)==34651
    assert next(r for r in targets if r['target']==67943)['known_positive_query_lower_bound']==1983744
    dest.mkdir();write_csv(dest/'effective_slot_plan_references.csv',ledger);write_csv(dest/'effective_target_cost.csv',targets);write_csv(dest/'remaining30_atlas_slots.csv',pending)
    helper=ROOT/'material_shared_collected_radius.h';inputs=[Path(__file__).resolve(),previous/'contract.json',previous/'effective_slot_plan_references.csv',coverage/'contract.json',run/'suffix_clearance.csv',path,policy,helper]
    (dest/'input.sha256').write_text(''.join(f'{sha(q)}  {q}\n' for q in inputs))
    result={'scope':'only adopted slot3200 unchanged actual atlas, saved checked prefix and same-function collected-ray tail proof; other failures preserved','newly_completed_atlas_slots':list(ready),'newly_adopted_panels':sum(r['panels'] for r in ready.values()),'active_new_panels':34651,'retained_original_target_panels':1047,'remaining_incomplete_atlas_slots':30,'authorized_target_query_limit':limit,'all_known_query_lower_bounds_within_new_limit':all(r['known_positive_query_lower_bound']<=limit for r in targets),'required_collected_radius_helper':str(helper),'required_collected_radius_helper_SHA256':sha(helper),'all_target_actual_atlas_complete':False,'all_target_analytic_certificates_complete':False,'full_quadrature_allowed':False,'new_physical_quadrature_queries':0,'target_cost':targets,'remaining_atlas_slots':pending}
    out=dest/'contract.json';out.write_text(json.dumps(result,indent=2)+'\n');(dest/'output.sha256').write_text(''.join(f'{sha(q)}  {q}\n' for q in (out,dest/'effective_slot_plan_references.csv',dest/'effective_target_cost.csv',dest/'remaining30_atlas_slots.csv')));print(json.dumps({k:result[k] for k in ('scope','newly_completed_atlas_slots','newly_adopted_panels','active_new_panels','remaining_incomplete_atlas_slots','authorized_target_query_limit','full_quadrature_allowed')},indent=2))
if __name__=='__main__':main()

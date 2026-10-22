"""Adopt exactly one independently covered/clearance-passed candidate plan."""
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
    dest=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_atlas_adoption_v2')
    assert dest.parent==ROOT and not dest.exists()
    previous=ROOT/'shared_atlas_adoption_v1';coverage=ROOT/'shared_origin_3208_coverage_v1';new=ROOT/'shared_origin_atlas_v2_mpi1'
    for run in (previous,coverage,new):
        for name in ('input.sha256','output.sha256'):
            for line in (run/name).read_text().splitlines():
                digest,path=line.split('  ',1);assert sha(Path(path))==digest
    covered=json.loads((coverage/'contract.json').read_text());assert covered['domain_coverage_and_native_geometry_density_clearance_complete'] and covered['slot_id']==3208
    policy=ROOT/'SHARED_ATLAS_BUDGET_POLICY_20261021.json';p=json.loads(policy.read_text());limit=p['new_per_target_query_limit'];assert limit==13074432
    data=new/'graded_source133_atlas_rank0.bin';raw=data.read_bytes();h=struct.unpack_from('<8Q',raw);assert h[3]==704
    indices=[j for j in range(h[4]) if struct.unpack_from('<Q',raw,64+j*704)[0]==3208];assert indices==list(range(3432,3704))
    ledger=rows(previous/'effective_slot_plan_references.csv');assert len(ledger)==9216
    costs=defaultdict(Counter);pending=[]
    for r in ledger:
        slot=int(r['slot_id'])
        if slot==3208:
            assert r['atlas_complete']=='0' and r['active_panels']=='0'
            r.update(atlas_complete='1',active_panels='272',complete_planned_queries='69632',provenance='adopted_source133_actual_graded_origin_plan',typed_packet=str(data),typed_record_first=str(indices[0]),typed_record_count='272')
        target=int(r['target']);v=costs[target];v['pending']+=r['atlas_complete']=='0'
        v['retained' if r['provenance']=='retained_original_target_plan' else 'new']+=int(r['active_panels'])
        if r['atlas_complete']=='0':pending.append({'slot_id':slot,'source':int(r['source']),'target':target,'original_failure':r['original_reason']})
    targets=[]
    for target,v in sorted(costs.items()):
        known=(v['new']+v['retained'])*256
        targets.append({'target':target,'incomplete_atlas_slots':v['pending'],'active_new_panels':v['new'],'retained_original_target_panels':v['retained'],'known_positive_query_lower_bound':known,'complete_planned_queries':'UNKNOWN' if v['pending'] else known,'authorized_query_cap':limit,'complete_cost_pass':int(not v['pending'] and known<=limit),'analytic_certificates_complete':0})
    assert len(pending)==32 and sum(r['active_new_panels'] for r in targets)==31251
    assert next(r for r in targets if r['target']==68527)['known_positive_query_lower_bound']==1360640
    dest.mkdir();write_csv(dest/'effective_slot_plan_references.csv',ledger);write_csv(dest/'effective_target_cost.csv',targets);write_csv(dest/'remaining32_atlas_slots.csv',pending)
    inputs=[Path(__file__).resolve(),previous/'contract.json',previous/'effective_slot_plan_references.csv',coverage/'contract.json',new/'origin_atlas_cost_and_clearance.csv',data,policy]
    (dest/'input.sha256').write_text(''.join(f'{sha(q)}  {q}\n' for q in inputs))
    contract={'scope':'add exactly covered and native-clearance-passed slot3208 plan to immutable effective registry under new user budget; failed candidates retained but not adopted','effective_slots':9216,'newly_completed_atlas_slot':3208,'newly_adopted_panels':272,'newly_planned_queries':69632,'active_new_panels':31251,'retained_original_target_panels':1047,'remaining_incomplete_atlas_slots':32,'authorized_target_query_limit':limit,'all_known_query_lower_bounds_within_new_budget':all(r['known_positive_query_lower_bound']<=limit for r in targets),'full_target_costs_unknown_count':sum(r['incomplete_atlas_slots']>0 for r in targets),'all_target_actual_atlas_complete':False,'all_target_analytic_certificates_complete':False,'new_physical_quadrature_queries':0,'full_quadrature_allowed':False,'target_cost':targets,'remaining_atlas_slots':pending}
    out=dest/'contract.json';out.write_text(json.dumps(contract,indent=2)+'\n');(dest/'output.sha256').write_text(''.join(f'{sha(q)}  {q}\n' for q in (out,dest/'effective_slot_plan_references.csv',dest/'effective_target_cost.csv',dest/'remaining32_atlas_slots.csv')))
    print(json.dumps({k:contract[k] for k in ('scope','newly_completed_atlas_slot','newly_adopted_panels','active_new_panels','remaining_incomplete_atlas_slots','authorized_target_query_limit','all_known_query_lower_bounds_within_new_budget','full_target_costs_unknown_count','full_quadrature_allowed')},indent=2))
if __name__=='__main__':main()

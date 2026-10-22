"""Expose fixed-plan query lower bounds even when full target costs are unknown."""
from collections import defaultdict
from pathlib import Path
import json
import sys
from audit_shared_regular_q_transfer_v1 import rows,sha
ROOT=Path(__file__).resolve().parent
def main():
    dest=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_atlas_cost_floor_v1')
    assert dest.parent==ROOT and not dest.exists()
    previous=ROOT/'shared_atlas_adoption_v1'
    for name in ('input.sha256','output.sha256'):
        for line in (previous/name).read_text().splitlines():
            digest,path=line.split('  ',1);assert sha(Path(path))==digest
    normal=ROOT/'SHARED_NORMALIZERS_20261021.json';n=json.loads(normal.read_text());assert n['other23_cold_normalizers_complete'] and n['saved_positive_subset_and_tau_accounting_pass']
    slots=rows(previous/'effective_slot_plan_references.csv');assert len(slots)==9216
    groups=defaultdict(list)
    for r in slots:groups[int(r['target'])].append(r)
    result=[]
    for target,records in sorted(groups.items()):
        complete=[r for r in records if r['atlas_complete']=='1']
        lower=sum(int(r['active_panels']) for r in complete)*256
        pending=sum(r['atlas_complete']=='0' for r in records)
        result.append({'target':target,'pending_atlas_slots':pending,'known_complete_plan_panels':lower//256,'known_positive_query_lower_bound':lower,'original_query_cap':817152,'already_exceeds_cap':lower>817152,'full_planned_queries':'UNKNOWN' if pending else lower,'fixed_existing_plans_can_meet_cap':False if lower>817152 else 'UNKNOWN' if pending else True,'largest_complete_contributors':[{'slot_id':int(r['slot_id']),'source':int(r['source']),'panels':int(r['active_panels'])} for r in sorted(complete,key=lambda r:-int(r['active_panels']))[:3]]})
    over=[r for r in result if r['already_exceeds_cap']]
    assert [(r['target'],r['known_positive_query_lower_bound']) for r in over]==[(67943,1182976),(68527,1291008)]
    assert [(r['largest_complete_contributors'][0]['source'],r['largest_complete_contributors'][0]['panels']) for r in over]==[(131,3078),(271,3078)]
    dest.mkdir();inputs=[Path(__file__).resolve(),previous/'contract.json',previous/'effective_slot_plan_references.csv',normal]
    (dest/'input.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in inputs))
    out={'scope':'exact immutable active complete-plan query lower bounds, no geometry/field/quad replay','already_over_cap_targets':over,'all_target_cost_floor':result,'two_complete_target_plan_families_need_replanning':True,'remaining_incomplete_atlas_slots':33,'all_target_atlas_complete':False,'full_quadrature_allowed':False,'new_physical_quadrature_queries':0}
    path=dest/'contract.json';path.write_text(json.dumps(out,indent=2)+'\n');(dest/'output.sha256').write_text(f'{sha(path)}  {path}\n')
    print(json.dumps({'scope':out['scope'],'already_over_cap_targets':over,'remaining_incomplete_atlas_slots':33,'full_quadrature_allowed':False},indent=2))
if __name__=='__main__':main()

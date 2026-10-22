"""Apply the new user-authorized limit to saved costs, preserving old verdicts."""
from pathlib import Path
import json
import sys
from audit_shared_regular_q_transfer_v1 import sha
ROOT=Path(__file__).resolve().parent
def main():
    dest=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_atlas_budget_v2')
    assert dest.parent==ROOT and not dest.exists()
    policy=ROOT/'SHARED_ATLAS_BUDGET_POLICY_20261021.json';p=json.loads(policy.read_text())
    assert p['new_per_target_query_limit']==16*817152==13074432 and p['new_per_target_panel_limit']*256==p['new_per_target_query_limit']
    old=ROOT/'shared_atlas_adoption_v1';floor=ROOT/'shared_atlas_cost_floor_v1';paths=[]
    for run in (old,floor):
        for name in ('input.sha256','output.sha256'):
            for line in (run/name).read_text().splitlines():
                digest,path=line.split('  ',1);assert sha(Path(path))==digest
        paths.append(run/'contract.json')
    source=json.loads((floor/'contract.json').read_text());targets=[];limit=p['new_per_target_query_limit']
    for r in source['all_target_cost_floor']:
        r=dict(r);known=r['known_positive_query_lower_bound'];pending=r['pending_atlas_slots']
        r.update({'previous_query_cap':r.pop('original_query_cap'),'new_query_cap':limit,'previous_known_lower_exceeded_cap':r.pop('already_exceeds_cap'),'known_lower_exceeds_new_cap':known>limit,'complete_cost_pass_under_new_cap':not pending and known<=limit,'full_cost_status_under_new_cap':'UNKNOWN' if pending else 'PASS' if known<=limit else 'EXCEEDS_NEW_CAP'})
        r.pop('fixed_existing_plans_can_meet_cap');targets.append(r)
    assert all(not r['known_lower_exceeds_new_cap'] for r in targets)
    assert sum(r['pending_atlas_slots'] for r in targets)==33
    review=ROOT/'ASTRA_SHARED_FAILED_ATLAS_BUDGET_ADDENDUM_20261021.md';assert review.is_file()
    paths.extend([Path(__file__).resolve(),policy,review])
    dest.mkdir();(dest/'input.sha256').write_text(''.join(f'{sha(q)}  {q}\n' for q in paths))
    contract={'scope':'new user-authorized16x query/panel limit applied ONLY to immutable saved plan references/costs','budget_policy_id':p['policy_id'],'budget_policy_SHA256':sha(policy),'budget_review_SHA256':sha(review),'new_per_target_query_limit':limit,'new_per_target_panel_limit':limit//256,'existing_completed_plans_preserved':True,'cost_only_replanning_of_sources131_271_needed':False,'known_lower_bounds_all_within_new_limit':True,'remaining_incomplete_atlas_slots':33,'full_target_costs_still_unknown_count':sum(r['pending_atlas_slots']>0 for r in targets),'target_cost':targets,'all_target_actual_atlas_complete':False,'all_target_analytic_certificates_complete':False,'full_quadrature_allowed':False,'new_physical_quadrature_queries':0,'streamed_all_target_quadrature_implemented':False}
    out=dest/'contract.json';out.write_text(json.dumps(contract,indent=2)+'\n');(dest/'output.sha256').write_text(f'{sha(out)}  {out}\n')
    print(json.dumps({k:contract[k] for k in ('scope','new_per_target_query_limit','new_per_target_panel_limit','cost_only_replanning_of_sources131_271_needed','known_lower_bounds_all_within_new_limit','remaining_incomplete_atlas_slots','full_target_costs_still_unknown_count','full_quadrature_allowed')},indent=2))
if __name__=='__main__':main()

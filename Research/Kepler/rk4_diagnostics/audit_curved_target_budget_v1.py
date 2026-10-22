"""Read-only reconstruction of positive normalizers and partial cold budgets.
No fluid evolution; the native interval/libm assumptions remain explicit.
"""
import csv, json, math, struct
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def power2(e):return F(1<<e) if e>=0 else F(1,1<<(-e))
def r80(s):
    x=F(s)
    if not x:return x
    sign=-1 if x<0 else 1;x=abs(x);e=x.numerator.bit_length()-x.denominator.bit_length()
    if x<power2(e):e-=1
    scaled=x*power2(63-e);n,rem=divmod(scaled.numerator,scaled.denominator)
    if rem*2>scaled.denominator or (rem*2==scaled.denominator and n%2):n+=1
    return sign*F(n)*power2(e-63)
def read(path):return list(csv.DictReader(open(ROOT/path)))
def dec(x):return D(x.numerator)/D(x.denominator)
def verify_regions(folder):
    rows=read(folder+'/interior_regions.csv');by={};ok=True
    for r in rows:
        target=int(r['target']);v=[(r80(r['x'+str(k)]),r80(r['y'+str(k)])) for k in range(3)]
        ux,uy=v[1][0]-v[0][0],v[1][1]-v[0][1];vx,vy=v[2][0]-v[0][0],v[2][1]-v[0][1]
        exact_area=(ux*vy-uy*vx)/2
        valid=(F(0)<r80(r['area_lower'])<=exact_area and r80(r['J_upper'])>0 and F(0)<=r80(r['lipschitz'])<1 and r80(r['owner_margin'])>0 and r80(r['U_lower'])>0 and r80(r['M_lower'])>0)
        box=[(r80(r['position_'+k+'lo']),r80(r['position_'+k+'hi'])) for k in ('x','y')]
        for old in by.get(target,[]):valid=valid and any(box[k][1]<old[k][0] or old[k][1]<box[k][0] for k in range(2))
        by.setdefault(target,[]).append(box);ok=ok and valid
    bounds={int(r['target']):r for r in read(folder+'/target_lower.csv')}
    for target,r in bounds.items():
        entries=[q for q in rows if int(q['target'])==target]
        ok=ok and r80(r['U_lower'])<=sum((r80(q['U_lower']) for q in entries),F(0)) and r80(r['M_lower'])<=sum((r80(q['M_lower']) for q in entries),F(0))
    return {'folder':folder,'regions':len(rows),'area_lower_and_image_disjointness_checks':ok,'native_status':(ROOT/folder/'status.txt').read_text().strip()},bounds

def main():
    report={'scope':'partial cold-transfer budget, not a full target certificate; bounds above budget do not prove physical error','normalizer_proofs':[],'cases':{},'missing_terms':['source_quadrature_and_source_reference_error','field_evaluation','analytic_trace_gluing','accumulation_rounding']}
    with localcontext() as ctx:
        ctx.prec=70
        for case in ('225','137','lens'):
            folders=['curved_target_lower_v1_'+case]+(['curved_target_lower_v2_lens'] if case=='lens' else [])
            chosen={}
            for folder in folders:
                check,bounds=verify_regions(folder);report['normalizer_proofs'].append(check)
                if not check['area_lower_and_image_disjointness_checks']:raise RuntimeError('proof data mismatch: '+folder)
                for target,r in bounds.items():
                    if D(r['U_lower'])>D(chosen.get(target,{}).get('U_lower','-1')):chosen[target]={**r,'proof_folder':folder}
            data=read('curved_fold_transfer_v5_'+case+'/target_contribution_bounds.csv');source=int(data[0]['source']);flow={}
            terms=['quad_estimate','representation_bound','native_boundary_bound','event_bound']
            for r in data:
                target=int(r['target'])
                if target==source:continue
                for affected in (source,target):
                    key=(affected,r['moment']);row=flow.setdefault(key,{k:D(0) for k in terms})
                    for k in terms:row[k]+=D(r[k])
            result=[]
            for (target,moment),values in sorted(flow.items()):
                if moment not in ('K','U','potential','thermal_rate'):continue
                proof=chosen[target];U=D(proof['U_lower']);share=D(proof['partial_budget_share']);tolerance=D('1e-9') if moment=='thermal_rate' else D('1e-10');budget=U*share*tolerance;total=sum(values.values(),D(0))
                result.append({'target':target,'moment':moment,'U_lower':str(U),'proof_folder':proof['proof_folder'],'fixed_partial_share':str(share),'partial_allowance':str(budget),'known_upper_bound_components':{k:str(v) for k,v in values.items()},'known_upper_bound_total':str(total),'known_bound_divided_by_U_lower':str(total/U),'known_bound_budget_ratio':str(total/budget),'remaining_terms_complete':False})
            report['cases'][case]={'source':source,'all_chosen_U_lower_positive':all(D(r['U_lower'])>0 for r in chosen.values()),'chosen_positive_proofs':chosen,'source_outgoing_and_receiver_incoming_share_same_record':True,'own_source_error_is_sum_of_outgoing_errors':True,'rows':result,'max_partial_known_budget_ratio':str(max((D(r['known_bound_budget_ratio']) for r in result),default=D(0))),'full_cold_transfer_proved':False}
            print(case,'max partial known budget ratio',report['cases'][case]['max_partial_known_budget_ratio'])
    path=ROOT/'CURVED_TARGET_PARTIAL_BUDGET_20261021.json'
    with open(path,'x') as f:json.dump(report,f,indent=2)
if __name__=='__main__':main()

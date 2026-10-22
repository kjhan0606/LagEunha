"""Exact preserved plans/prefixes and collected-ray suffix proof accounting."""
from collections import Counter,defaultdict
from fractions import Fraction as F
from pathlib import Path
import json
import struct
import sys
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_regular_q_transfer_v1 import rows,sha
ROOT=Path(__file__).resolve().parent
def ball(b,at):return quad(b,at),quad(b,at+16)
def main():
    run=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_ray_suffix_v1_mpi1');assert run.parent==ROOT
    for name in ('source.sha256','input.sha256','binary.sha256','output.sha256'):
        for line in (run/name).read_text().splitlines():
            digest,path=line.split('  ',1);assert sha(Path(path))==digest
    for name in ('retained_source_integrity.txt','retained_output_integrity.txt','final_source_integrity.txt','final_input_integrity.txt'):
        lines=(run/name).read_text().splitlines();assert lines and all(s.endswith(': OK') for s in lines)
    assert (run/'native_exit_code.txt').read_text().strip()=='0';c=json.loads((run/'contract.json').read_text())
    oldrun=ROOT/'shared_origin_atlas_v2_mpi1';_,old=wire(oldrun/'graded_source133_atlas_rank0.bin')
    oldrows={(struct.unpack_from('<Q',b)[0],struct.unpack_from('<Q',b,16)[0]):b for b in old if struct.unpack_from('<Q',b)[0] in (3200,3203)}
    _,plans=wire(run/'unchanged_original_candidate_plans_rank0.bin');assert len(plans)==3432 and all(b==oldrows[struct.unpack_from('<Q',b)[0],struct.unpack_from('<Q',b,16)[0]] for b in plans)
    byplan={(struct.unpack_from('<Q',b)[0],struct.unpack_from('<Q',b,16)[0]):b for b in plans}
    first={3200:94,3203:2};total={3200:3128,3203:304}
    ledger={int(r['slot_id']):r for r in rows(run/'suffix_clearance.csv')};assert set(ledger)==set(first)
    oldledger={int(r['slot_id']):r for r in rows(oldrun/'origin_atlas_cost_and_clearance.csv')}
    h,proof=wire(run/'suffix_axis_proofs_rank0.bin');assert h[3]==592
    byslot=defaultdict(list);poly_checks=0
    for b in proof:
        slot,source,target,index,axis,passed=struct.unpack_from('<6Q',b);assert source==133 and first[slot]<=index<total[slot] and axis in (0,1)
        byslot[slot].append((index,axis,passed));assert target==int(ledger[slot]['target'])
        p=byplan[slot,index];chart=struct.unpack_from('<i',p,48)[0]
        if passed:
            contraction,invariance=quad(b,304),quad(b,320);assert 0<=contraction<1 and 0<=invariance<1
            if chart in (2,3):
                n,e=ball(b,144);argument=ball(b,240);assert argument[0]>argument[1]>=0
                if chart==2:center=(9+n*n)/1048576;radius=(2*abs(n)*e+e*e)/1048576
                else:center=(40-32*n+10*n*n)/4194304;radius=(abs(-32+20*n)*e+10*e*e)/4194304
                assert argument[0]-argument[1]<=center-radius<=center+radius<=argument[0]+argument[1]
                poly_checks+=1
    completed=0
    for slot,r in ledger.items():
        assert int(r['total_original_panels'])==total[slot] and int(r['reused_prefix_axis_boxes'])==int(oldledger[slot]['checked_axis_boxes'])==first[slot]*2
        v=byslot[slot];assert len(v)==int(r['new_assessed_axis_boxes'])
        assert [(i,a) for i,a,_ in v]==[(first[slot]+j//2,j%2) for j in range(len(v))]
        if r['clearance_complete']=='1':
            completed+=1;assert len(v)==2*(total[slot]-first[slot]) and all(passed for _,_,passed in v)
            assert r['first_failed_panel']==r['first_failed_axis']=='-1'
        else:
            failures=[(i,a) for i,a,p in v if not p];assert failures and min(failures)==(int(r['first_failed_panel']),int(r['first_failed_axis']))
            assert len(v)%32==0 or first[slot]+len(v)//2==total[slot]
        assert r['query_cap']=='13074432'
    assert completed==c['completed_candidate_geometry_density_clearance'] and len(proof)==c['new_assessed_axis_boxes'] and c['reused_prefix_axis_boxes']==192
    assert c['new_physical_quadrature_queries']==c['new_inverse_solves']==c['new_root_jobs']==c['normal_panel_replays']==0 and c['state_unchanged']
    assert not c['all12_analytic_moment_remainders_complete'] and not c['full_quadrature_allowed']
    print(json.dumps({'scope':'exact unchanged original plans/passed-prefix accounting and new suffix/collected-ray polynomial disk inequalities; retained implicit graph/field engineering premises','run':str(run),'contract':c,'independent_preserved_plan_and_suffix_accounting_pass':True,'independent_collected_ray_polynomial_disk_checks':poly_checks,'candidate_domain_coverage_not_yet_independently_adopted':True,'slot_results':list(ledger.values()),'new_physical_quadrature_queries':0,'full_quadrature_allowed':False},indent=2))
if __name__=='__main__':main()

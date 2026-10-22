"""Saved graded candidate counts/provenance/first-clearance-failure audit."""
from collections import Counter
from pathlib import Path
import json
import struct
import sys
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_regular_q_transfer_v1 import Q,rows,sha
ROOT=Path(__file__).resolve().parent
def ball(b,at):return quad(b,at),quad(b,at+16)
def main():
    run=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_origin_atlas_v2_mpi1');assert run.parent==ROOT
    for name in ('source.sha256','input.sha256','binary.sha256','output.sha256'):
        for line in (run/name).read_text().splitlines():
            digest,path=line.split('  ',1);assert sha(Path(path))==digest
    for name in ('retained_source_integrity.txt','retained_arrangement_integrity.txt','final_source_integrity.txt','final_input_integrity.txt'):
        lines=(run/name).read_text().splitlines();assert lines and all(s.endswith(': OK') for s in lines)
    assert (run/'native_exit_code.txt').read_text().strip()=='0'
    c=json.loads((run/'contract.json').read_text());policy=json.loads((ROOT/'SHARED_ATLAS_BUDGET_POLICY_20261021.json').read_text())
    assert c['authorized_target_query_limit']==policy['new_per_target_query_limit']==13074432
    ledger=rows(run/'origin_atlas_cost_and_clearance.csv');expected_ids={int(x) for x in sys.argv[2].split(',')} if len(sys.argv)>2 else {3200,3203,3208,3211};assert {int(r['slot_id']) for r in ledger}==expected_ids
    byslot={int(r['slot_id']):r for r in ledger};counts=Counter();whole=Counter();clipped=Counter()
    _,panels=wire(run/'graded_source133_atlas_rank0.bin')
    for b in panels:
        slot,source,index,target,policy,dep,lower,upper=struct.unpack_from('<4Q4i',b)
        assert source==133 and slot in byslot and target==int(byslot[slot]['target']) and index==counts[slot];counts[slot]+=1
        ch,axis,lo,hi=struct.unpack_from('<4i',b,48);assert 0<=ch<6 and (lo,hi)==(lower,upper)
        assert quad(b,672)>quad(b,688)>=0
        iswhole=b[192]!=0
        assert policy==(0 if iswhole else 1 if dep==1 else 2)
        if iswhole:
            whole[slot]+=1;chart,ir,it,nr,nt=struct.unpack_from('<5i',b,208)
            assert chart==ch and 0<=ir<nr and 0<=it<nt and nt>=4 and dep==-1
            for a,z in ((240,272),(304,336)):assert quad(b,a)<quad(b,z)
        else:
            clipped[slot]+=1;assert dep==axis and dep in (0,1) and -1<=lo<23 and -1<=hi<23
            assert quad(b,64)<quad(b,96) and quad(b,128)<quad(b,160)
        for at in ((256,288,320,352) if iswhole else (80,112,144,176)):assert quad(b,at)>=0
    planned=clear=checked=0
    for r in ledger:
        slot=int(r['slot_id']);n=int(r['stored_panels']);assert n==counts[slot] and n*256<=13074432
        p=r['plan_complete']=='1';ok=r['clearance_complete']=='1';planned+=p;clear+=ok;axisboxes=int(r['checked_axis_boxes']);checked+=axisboxes
        assert r['predicted_queries']==(str(n*256) if p else 'UNKNOWN')
        assert r['authorized_query_cap']=='13074432'
        if ok:assert p and axisboxes==2*n and r['failed_panel']==r['failed_axis']=='-1'
        elif p:
            i=int(r['failed_panel']);axis=int(r['failed_axis']);assert 0<=i<n and axis in (0,1) and axisboxes==2*i+axis
        else:assert axisboxes==0 and r['failed_panel']==r['failed_axis']=='-1'
        assert r['reason']
    assert len(panels)==c['stored_panels'] and planned==c['completed_fixed_plans'] and clear==c['completed_geometry_density_clearance'] and checked==c['checked_axis_boxes']
    assert c['retained_normal_chart_certificates']==16 and c['normal_chart_reconstructions']==c['new_root_jobs']==c['new_inverse_solves']==c['new_physical_quadrature_queries']==0
    assert c['state_unchanged'] and not c['all_target_analytic_certificates_complete'] and not c['full_quadrature_allowed']
    print(json.dumps({'scope':'saved graded candidate ABI/identity/count/first-failure accounting only; not an independent graded coverage or complex field proof','run':str(run),'contract':c,'independent_candidate_and_failure_accounting_pass':True,'slot_witnesses':ledger,'ordinary_panel_counts':dict(whole),'implicit_panel_counts':dict(clipped),'graded_actual_domain_coverage_independently_adopted':False,'all_target_analytic_certificates_complete':False,'full_quadrature_allowed':False},indent=2))
if __name__=='__main__':main()

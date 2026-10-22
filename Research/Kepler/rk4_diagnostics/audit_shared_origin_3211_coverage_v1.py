"""Exact partition proof for the one focused radial-endpoint correction."""
from collections import defaultdict
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
    dest=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_origin_3211_coverage_v1');assert dest.parent==ROOT and not dest.exists()
    audit=ROOT/'SHARED_ORIGIN_ATLAS_V3_20261021.json';a=json.loads(audit.read_text());assert a['independent_candidate_and_failure_accounting_pass']
    events=ROOT/'shared_inverse_events_v1_mpi1';run=ROOT/'shared_origin_atlas_v3_mpi1'
    _,v=wire(events/'new_partial_strips_rank0.bin');base=[]
    for b in v:
        slot,source,target,chart,index=struct.unpack_from('<5Q',b)
        if slot!=3211:continue
        ch,dep,lower,upper=struct.unpack_from('<4i',b,48);assert source==133 and target==68584 and ch==chart
        base.append({'key':(chart,index),'chart':chart,'dep':dep,'lower':lower,'upper':upper,'lo':ball(b,64),'hi':ball(b,96)})
    assert len(base)==5 and sum(r['dep']==0 for r in base)==1
    old=ROOT/'shared_all_atlas_v2_mpi4';excluded=[r for rank in range(4) for r in rows(old/f'actual_chart_status_rank{rank}.csv') if r['slot_id']=='3211' and r['chart'] not in ('2','3')];assert len(excluded)==4 and all(r['complete']==r['empty']=='1' for r in excluded)
    _,v=wire(events/'new_failed_chart_arrangements_rank0.bin');reference={}
    for b in v:
        slot,source,target,chart,complete=struct.unpack_from('<5Q',b)
        if slot==3211:assert complete;reference[chart]=b[96:432]
    r=[r for r in rows(run/'origin_atlas_cost_and_clearance.csv') if r['slot_id']=='3211'];assert len(r)==1;r=r[0]
    assert r['plan_complete']==r['clearance_complete']=='1' and r['checked_axis_boxes']=='544'
    _,v=wire(run/'graded_source133_atlas_rank0.bin');assert len(v)==272
    groups=defaultdict(lambda:defaultdict(list));policy_counts=defaultdict(int)
    for index,b in enumerate(v):
        slot,source,j,target,policy,dep,lower,upper=struct.unpack_from('<4Q4i',b);assert (slot,source,j,target)==(3211,133,index,68584)
        ch,axis,lo,hi=struct.unpack_from('<4i',b,48);assert b[368:704]==reference[ch]
        if b[192]:
            assert policy==0 and dep==-1 and lo==hi==-1
            free_left,free_right=ball(b,240),ball(b,272);eta_left,eta_right=ball(b,304),ball(b,336);saved_dep=1;lower=upper=-1
        else:
            assert dep==axis and policy==(1 if dep else 2)
            free_left,free_right=ball(b,64),ball(b,96);eta_left,eta_right=ball(b,128),ball(b,160);saved_dep=dep
        candidates=[p for p in base if (p['chart'],p['dep'],p['lower'],p['upper'])==(ch,saved_dep,lower,upper) and p['lo'][0]<=free_left[0]<free_right[0]<=p['hi'][0]]
        assert len(candidates)==1;parent=candidates[0]
        groups[parent['key']][free_left,free_right].append((eta_left,eta_right));policy_counts[policy]+=1
    free_parts=eta_parts=0
    for p in base:
        children=sorted(groups[p['key']].items(),key=lambda r:r[0][0][0]);assert children
        if p['dep']==1:
            assert children[0][0][0]==p['lo'] and children[-1][0][1]==p['hi']
        else:
            # The recorded angular endpoints enclose the declared exact
            # affine j/16 values. Ball evaluation of 0+1*0 adds outward
            # roundoff; it does not change that mathematical endpoint.
            assert p['lo']==(0,0) and p['hi']==(1,0) and len(children)==16
            for j,((lo,hi),_) in enumerate(children):
                assert lo[0]-lo[1]<=F(j,16)<=lo[0]+lo[1]
                assert hi[0]-hi[1]<=F(j+1,16)<=hi[0]+hi[1]
        for ((_,hi),_),((lo,_),_) in zip(children,children[1:]):assert hi==lo
        free_parts+=len(children)
        for (lo,hi),etas in children:
            assert hi[0]-lo[0]>hi[1]+lo[1]
            etas.sort(key=lambda r:r[0][0]);eta_parts+=len(etas)
            assert etas[0][0][0]==0 and etas[-1][1][0]==1
            for (lo,hi) in etas:assert hi[0]-lo[0]>hi[1]+lo[1]
            for (_,hi),(lo,_) in zip(etas,etas[1:]):assert hi==lo
            if p['dep']==1:
                n=len(etas);assert n>=16
                for j,(lo,hi) in enumerate(etas):
                    assert lo[0]-lo[1]<=F(j,n)<=lo[0]+lo[1]
                    assert hi[0]-hi[1]<=F(j+1,n)<=hi[0]+hi[1]
            else:
                assert all(lo[1]==hi[1]==0 for lo,hi in etas)
                # Positive exact Q128 differences include the small final
                # interval near1. No endpoint is snapped or discarded.
    assert eta_parts==272
    dest.mkdir();paths=[Path(__file__).resolve(),audit,events/'new_partial_strips_rank0.bin',events/'new_failed_chart_arrangements_rank0.bin',run/'graded_source133_atlas_rank0.bin',run/'origin_atlas_cost_and_clearance.csv',*(old/f'actual_chart_status_rank{k}.csv' for k in range(4))]
    (dest/'input.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in paths))
    result={'scope':'exact saved source133/target68584 original-strip partition and native-clearance accounting, includes exact tiny eta endpoint interval; no independent field or analytic moment proof','slot_id':3211,'source':133,'target':68584,'actual_original_strips':5,'graded_free_subslabs':free_parts,'panels':272,'queries':69632,'panel_policies':dict(policy_counts),'all_original_event_endpoints_and_axes_preserved':True,'tiny_endpoint_interval_retained_not_snapped':True,'domain_coverage_and_native_geometry_density_clearance_complete':True,'all12analytic_moment_remainders_complete':False,'new_physical_quadrature_queries':0,'full_quadrature_allowed':False}
    path=dest/'contract.json';path.write_text(json.dumps(result,indent=2)+'\n');(dest/'output.sha256').write_text(f'{sha(path)}  {path}\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()

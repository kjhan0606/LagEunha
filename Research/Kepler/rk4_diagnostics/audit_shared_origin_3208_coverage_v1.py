"""Exact saved interval/branch partition audit of the passing slot3208 plan."""
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
    dest=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_origin_3208_coverage_v1')
    assert dest.parent==ROOT and not dest.exists()
    audit=ROOT/'SHARED_ORIGIN_ATLAS_20261021.json';a=json.loads(audit.read_text());assert a['independent_candidate_and_failure_accounting_pass']
    events=ROOT/'shared_inverse_events_v1_mpi1';run=ROOT/'shared_origin_atlas_v2_mpi1'
    _,strips=wire(events/'new_partial_strips_rank0.bin');base=[]
    for b in strips:
        slot,source,target,chart,index=struct.unpack_from('<5Q',b)
        if slot!=3208:continue
        ch,dep,lower,upper=struct.unpack_from('<4i',b,48);assert source==133 and target==68527 and chart==ch==2 and dep==1
        base.append({'index':index,'chart':chart,'dep':dep,'lower':lower,'upper':upper,'lo':ball(b,64),'hi':ball(b,96)})
    assert len(base)==4 and [r['index'] for r in base]==list(range(4))
    _,charts=wire(events/'new_failed_chart_arrangements_rank0.bin');reference={}
    for b in charts:
        slot,source,target,chart,complete=struct.unpack_from('<5Q',b)
        if slot==3208:assert complete;reference[chart]=b[96:432]
    old=ROOT/'shared_all_atlas_v2_mpi4';oldstatuses=[r for rank in range(4) for r in rows(old/f'actual_chart_status_rank{rank}.csv') if r['slot_id']=='3208' and r['chart'] not in ('2','3')]
    assert len(oldstatuses)==4 and all(r['complete']==r['empty']=='1' for r in oldstatuses)
    # Replacement chart3 also contains no strips. All actual positive pieces
    # therefore are precisely the four adopted original-coordinate strips.
    chartrows=[r for r in rows(events/'inverse_chart_arrangements.csv') if r['slot_id']=='3208'];assert len(chartrows)==2
    assert {int(r['chart']):int(r['strips']) for r in chartrows}=={2:4,3:0}
    ledger=[r for r in rows(run/'origin_atlas_cost_and_clearance.csv') if r['slot_id']=='3208'];assert len(ledger)==1;ledger=ledger[0]
    assert ledger['plan_complete']==ledger['clearance_complete']=='1' and ledger['checked_axis_boxes']=='544'
    _,panels=wire(run/'graded_source133_atlas_rank0.bin');groups=defaultdict(lambda:defaultdict(list));count=0
    for b in panels:
        slot,source,index,target,policy,dep,lower,upper=struct.unpack_from('<4Q4i',b)
        if slot!=3208:continue
        assert source==133 and target==68527 and index==count;count+=1
        assert policy==1 and dep==1 and b[192]==0
        chart,axis,lo,hi=struct.unpack_from('<4i',b,48);assert chart==2 and axis==dep and (lo,hi)==(lower,upper)
        assert b[368:704]==reference[chart]
        left,right=ball(b,64),ball(b,96)
        matches=[r for r in base if (r['lower'],r['upper'])==(lower,upper) and r['lo'][0]<=left[0]<right[0]<=r['hi'][0]]
        assert len(matches)==1
        groups[matches[0]['index']][left,right].append((ball(b,128),ball(b,160)))
    assert count==272
    child_count=0
    for r in base:
        children=sorted(groups[r['index']].items(),key=lambda pair:pair[0][0][0]);child_count+=len(children)
        assert children and children[0][0][0]==r['lo'] and children[-1][0][1]==r['hi']
        for ((_,right),_),((left,_),_) in zip(children,children[1:]):assert right==left
        for (left,right),etas in children:
            assert left[0]+left[1]<right[0]-right[1]
            etas.sort(key=lambda pair:pair[0][0]);n=len(etas);assert n>=16
            for j,(lo,hi) in enumerate(etas):
                assert lo[0]-lo[1]<=F(j,n)<=lo[0]+lo[1]
                assert hi[0]-hi[1]<=F(j+1,n)<=hi[0]+hi[1]
            for (_,hi),(lo,_) in zip(etas,etas[1:]):assert hi==lo
    assert child_count==17
    dest.mkdir();paths=[Path(__file__).resolve(),audit,events/'new_failed_chart_arrangements_rank0.bin',events/'new_partial_strips_rank0.bin',events/'inverse_chart_arrangements.csv',run/'graded_source133_atlas_rank0.bin',run/'origin_atlas_cost_and_clearance.csv',*(old/f'actual_chart_status_rank{k}.csv' for k in range(4))]
    (dest/'input.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in paths))
    result={'scope':'exact branch/event-slab/eta partition of saved slot3208 candidate under adopted complete original-coordinate arrangement; native complex/pole clearance engineering premises retained','slot_id':3208,'source':133,'target':68527,'actual_original_strips':4,'graded_radial_child_slabs':17,'panels':272,'queries':69632,'all_free_boundaries_preserve_exact_saved_event_identity':True,'every_eta_partition_preserves_original_domain':True,'ordinary_chart_exclusions_reused':4,'original_chart3_zero_strip_arrangement_reused':True,'domain_coverage_and_native_geometry_density_clearance_complete':True,'all12analytic_moment_remainders_complete':False,'new_physical_quadrature_queries':0,'full_quadrature_allowed':False}
    out=dest/'contract.json';out.write_text(json.dumps(result,indent=2)+'\n');(dest/'output.sha256').write_text(f'{sha(out)}  {out}\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()

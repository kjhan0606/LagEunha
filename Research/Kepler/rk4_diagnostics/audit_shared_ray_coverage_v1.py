"""Exact six-chart/strip partition for completed same-ray candidate tails."""
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
def contain(b,x):return b[0]-b[1]<=x<=b[0]+b[1]
def main():
    dest=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_ray_coverage_v1');assert dest.parent==ROOT and not dest.exists()
    audit=ROOT/'SHARED_RAY_SUFFIX_20261021.json';a=json.loads(audit.read_text());assert a['independent_preserved_plan_and_suffix_accounting_pass']
    run=ROOT/'shared_ray_suffix_v1_mpi1';events=ROOT/'shared_inverse_events_v1_mpi1';old=ROOT/'shared_all_atlas_v2_mpi4'
    completed={int(r['slot_id']) for r in rows(run/'suffix_clearance.csv') if r['clearance_complete']=='1'};assert completed
    _,raw=wire(events/'new_partial_strips_rank0.bin');base=defaultdict(list)
    for b in raw:
        slot,source,target,chart,index=struct.unpack_from('<5Q',b)
        if slot not in completed:continue
        ch,dep,lower,upper=struct.unpack_from('<4i',b,48);assert ch==chart
        base[slot].append({'key':(chart,index),'chart':chart,'dep':dep,'lower':lower,'upper':upper,'lo':ball(b,64),'hi':ball(b,96)})
    _,raw=wire(events/'new_failed_chart_arrangements_rank0.bin');reference={}
    for b in raw:
        slot,source,target,chart,complete=struct.unpack_from('<5Q',b)
        if slot in completed:assert complete;reference[slot,chart]=b[96:432]
    statuses={(int(r['slot_id']),int(r['chart'])):r for rank in range(4) for r in rows(old/f'actual_chart_status_rank{rank}.csv') if int(r['slot_id']) in completed}
    _,panels=wire(run/'unchanged_original_candidate_plans_rank0.bin');by=defaultdict(list)
    for b in panels:
        slot=struct.unpack_from('<Q',b)[0]
        if slot in completed:by[slot].append(b)
    results=[]
    for slot,records in sorted(by.items()):
        children=defaultdict(lambda:defaultdict(list));normal=defaultdict(set);normal_sizes={};ordinary=implicit=0
        for index,b in enumerate(records):
            sid,source,j,target,policy,dep,lower,upper=struct.unpack_from('<4Q4i',b);assert sid==slot and source==133 and j==index
            chart,axis,lo,hi=struct.unpack_from('<4i',b,48)
            if chart not in (2,3):
                assert statuses[slot,chart]['complete']==statuses[slot,chart]['whole']=='1' and b[192] and policy==0
                ch,ir,it,nr,nt=struct.unpack_from('<5i',b,208);assert chart==ch
                if chart in normal_sizes:assert normal_sizes[chart]==(nr,nt)
                normal_sizes[chart]=nr,nt;assert (ir,it) not in normal[chart];normal[chart].add((ir,it));ordinary+=1
                for at,v in zip((240,272,304,336),(F(ir,nr),F(ir+1,nr),F(it,nt),F(it+1,nt))):assert contain(ball(b,at),v)
                continue
            assert b[368:704]==reference[slot,chart]
            if b[192]:
                assert policy==0 and dep==-1 and lo==hi==-1;left,right=ball(b,240),ball(b,272);eta0,eta1=ball(b,304),ball(b,336);saved_dep=1;lower=upper=-1
            else:
                assert dep==axis and policy==(1 if dep else 2);left,right=ball(b,64),ball(b,96);eta0,eta1=ball(b,128),ball(b,160);saved_dep=dep;implicit+=1
            possible=[p for p in base[slot] if (p['chart'],p['dep'],p['lower'],p['upper'])==(chart,saved_dep,lower,upper) and p['lo'][0]<=left[0]<right[0]<=p['hi'][0]]
            assert len(possible)==1;children[possible[0]['key']][left,right].append((eta0,eta1))
        for chart in (0,1,4,5):
            r=statuses[slot,chart];assert r['complete']=='1' and r['strips']==r['root_windows']=='0'
            if r['empty']=='1':assert not normal[chart]
            else:
                nr,nt=normal_sizes[chart];assert normal[chart]=={(i,j) for i in range(nr) for j in range(nt)}
        free_count=0;eta_count=0
        for p in base[slot]:
            v=sorted(children[p['key']].items(),key=lambda r:r[0][0][0]);assert v;free_count+=len(v)
            if p['dep']==1:assert v[0][0][0]==p['lo'] and v[-1][0][1]==p['hi']
            else:
                assert p['lo']==(0,0) and p['hi']==(1,0) and len(v)==16
                for j,((lo,hi),_) in enumerate(v):assert contain(lo,F(j,16)) and contain(hi,F(j+1,16))
            for ((_,hi),_),((lo,_),_) in zip(v,v[1:]):assert hi==lo
            for (lo,hi),etas in v:
                assert hi[0]-lo[0]>hi[1]+lo[1];etas.sort(key=lambda r:r[0][0]);eta_count+=len(etas)
                assert etas[0][0][0]==0 and etas[-1][1][0]==1
                for lo,hi in etas:assert hi[0]-lo[0]>hi[1]+lo[1]
                for (_,hi),(lo,_) in zip(etas,etas[1:]):assert hi==lo
                if p['dep']==1:
                    n=len(etas);assert n>=16
                    for j,(lo,hi) in enumerate(etas):assert contain(lo,F(j,n)) and contain(hi,F(j+1,n))
                else:assert all(lo[1]==hi[1]==0 for lo,hi in etas)
        assert eta_count+ordinary==len(records)
        results.append({'slot_id':slot,'source':133,'target':struct.unpack_from('<Q',records[0],24)[0],'panels':len(records),'queries':len(records)*256,'original_actual_strips':len(base[slot]),'graded_free_subslabs':free_count,'normal_whole_grid_panels':ordinary,'implicit_panels':implicit,'exact_same_event_branch_domain_partition_pass':True,'saved_prefix_plus_new_collected_ray_clearance_complete':True,'all12_analytic_moment_remainders_complete':False})
    dest.mkdir();paths=[Path(__file__).resolve(),audit,run/'suffix_clearance.csv',run/'unchanged_original_candidate_plans_rank0.bin',events/'new_partial_strips_rank0.bin',events/'new_failed_chart_arrangements_rank0.bin',*(old/f'actual_chart_status_rank{k}.csv' for k in range(4))]
    (dest/'input.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in paths));result={'scope':'exact unchanged candidate six-chart/strip/eta domain partition; retained arrangement and native complex-field engineering premises','adoptable_slot_plans':results,'domain_coverage_and_native_geometry_density_clearance_complete':True,'all12_analytic_moment_remainders_complete':False,'new_physical_quadrature_queries':0,'full_quadrature_allowed':False}
    out=dest/'contract.json';out.write_text(json.dumps(result,indent=2)+'\n');(dest/'output.sha256').write_text(f'{sha(out)}  {out}\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()

"""Exact pending29 rectangle coverage, geometric distances and full24 cost.

Audits a fixed candidate plan, not its complex clearance or all12 remainder.
"""
from collections import Counter,defaultdict
from fractions import Fraction as F
from pathlib import Path
import json,struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball,add,mul,neg,encloses,verify
from audit_shared_regular_q_transfer_v1 import rows,sha

ROOT=Path(__file__).resolve().parent

def main():
    run=ROOT/'shared_pending_atlas_v1_mpi1';c=verify(run)
    for name in ('registry_input_integrity.txt','registry_output_integrity.txt'):
        v=(run/name).read_text().splitlines();assert v and all(s.endswith(': OK') for s in v)
    ledger={int(r['slot_id']):r for r in rows(run/'pending29_typed_plan_cost.csv')};assert len(ledger)==29
    binding=rows(ROOT/'shared_pending29_chart_assembly_v1/actual174_chart_bindings.csv')
    cache={};references={};kind={};parents=[]
    _,strips=wire(ROOT/'shared_missing_arrangements_v1_mpi1/recovered_actual_strips_rank0.bin')
    for r in binding:
        key=int(r['slot_id']),int(r['chart']);kind[key]=r['kind']
        if r['kind']=='oblique':continue
        path=Path(r['typed_packet'])
        if path not in cache:cache[path]=wire(path)[1]
        b=cache[path][int(r['record'])];references[key]=b
        if r['kind']=='whole':parents.append(dict(key=(*key,-1),chart=key,dep=1,lower=-1,upper=-1,lo=(F(0),F(0)),hi=(F(1),F(0)),synthetic=True))
    for b in strips:
        slot,source,target,ch,index=struct.unpack_from('<5Q',b);ci,dep,lower,upper=struct.unpack_from('<4i',b,48);assert ci==ch
        parents.append(dict(key=(slot,ch,index),chart=(slot,ch),dep=dep,lower=lower,upper=upper,lo=ball(b,64),hi=ball(b,96),synthetic=False))
    assert len(parents)==205
    bychart=defaultdict(list)
    for p in parents:bychart[p['chart']].append(p)
    hm,metrics=wire(run/'actual_chart_density_metrics_rank0.bin');assert hm[3]==352 and len(metrics)==124
    method={};state=(ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin').read_bytes()
    with (ROOT/'material_p3_finite_v2/mpi1_cpu/curved_map_control.bin').open('rb') as f:f.read(64);parameters=[F(x) for x in struct.unpack('<12d',f.read(96))]
    for b in metrics:
        slot,source,target,ch,mode=struct.unpack_from('<5Q',b);key=slot,ch;z=references[key];assert kind[key]!='excluded' and key not in method;method[key]=mode
        assert (source,target)==struct.unpack_from('<2Q',z,8)
        v=[tuple(F(x) for x in struct.unpack_from('<2d',state,64+1064*source+16*j)) for j in range(3)]
        apex,part=struct.unpack_from('<2i',z,96)
        weights=[[F(int(j==apex)) for j in range(3)],
                 [F(1,3)]*3 if part else [F(int(j in (apex,(apex+1)%3)),2) for j in range(3)],
                 [F(int(j in (apex,(apex+2)%3)),2) for j in range(3)] if part else [F(1,3)]*3]
        P=[tuple(sum(weights[a][j]*v[j][k] for j in range(3)) for k in range(2)) for a in range(3)]
        lo,hi,scale=quad(b,304),quad(b,320),quad(b,336);assert lo>0
        if mode==0:
            for a in range(3):
                for k in range(2):encloses(ball(b,48+64*a+32*k),(P[a][k],F(0)))
            N=[ball(b,240+32*k)[0] for k in range(2)]
            sep=min(sum(N[k]*ball(b,48+64*a+32*k)[0]-abs(N[k])*ball(b,48+64*a+32*k)[1] for k in range(2)) for a in range(3))
            square=sum(x*x for x in N)
            assert square>0 and sep>0 and hi*hi>=square and lo*lo*square<=sep*sep
            assert 0<scale<=min(lo,parameters[3])/8
        else:
            assert mode==1 and P[0]==(0,0) and scale==parameters[7]
            U=P[1];V=tuple(P[2][k]-P[1][k] for k in range(2))
            a=sum(x*x for x in U);bb=2*sum(U[k]*V[k] for k in range(2));cc=sum(x*x for x in V)
            candidates=[a,a+bb+cc]
            if cc and 0<-bb/(2*cc)<1:candidates.append(a-bb*bb/(4*cc))
            assert lo*lo<=min(candidates) and hi*hi>=max(a,a+bb+cc)
    hp,panels=wire(run/'new_regular_typed_panels_rank0.bin');assert hp[3]==704 and len(panels)==50782
    _,oblique=wire(run/'unchanged_oblique_panels_rank0.bin');assert oblique==wire(ROOT/'shared_oblique_collected_suffix_v1_mpi1/unchanged_candidate_panels_rank0.bin')[1]
    old_oblique=Counter(struct.unpack_from('<Q',b)[0] for b in oblique)
    groups=defaultdict(lambda:defaultdict(list));counts=Counter();policy=Counter()
    for b in panels:
        slot,source,index,target,pol,dep,lower,upper=struct.unpack_from('<4Q4i',b);ch,axis,l,u=struct.unpack_from('<4i',b,48)
        assert slot in ledger and source==int(ledger[slot]['source']) and target==int(ledger[slot]['target']) and index==counts[slot]
        counts[slot]+=1;policy[pol]+=1;key=slot,ch;assert key in references and kind[key] not in ('oblique','excluded') and b[368:704]==references[key][96:432]
        if b[192]:
            assert pol==0 and dep==-1 and l==u==lower==upper==-1
            fl,fr,el,er=ball(b,240),ball(b,272),ball(b,304),ball(b,336);saved_dep=1
            ci,ir,it,nr,nt=struct.unpack_from('<5i',b,208);assert ci==ch and 0<=ir<nr and 0<=it<nt
        else:
            assert dep==axis and pol==(1 if dep else 2) and (l,u)==(lower,upper)
            fl,fr,el,er=ball(b,64),ball(b,96),ball(b,128),ball(b,160);saved_dep=dep
        assert fl[0]+fl[1]<fr[0]-fr[1] and el[0]+el[1]<er[0]-er[1]
        candidates=[p for p in bychart[key] if (p['dep'],p['lower'],p['upper'])==(saved_dep,lower,upper) and p['lo'][0]<=fl[0]<fr[0]<=p['hi'][0]]
        assert len(candidates)==1,(slot,index,ch,len(candidates));parent=candidates[0]
        groups[parent['key']][fl,fr].append((el,er))
    free_children=0
    for p in parents:
        children=sorted(groups[p['key']].items(),key=lambda r:r[0][0][0]);assert children,p['key']
        if method[p['chart']]==1 and p['dep']==0:
            assert len(children)==16
            width=add(p['hi'],neg(p['lo']))
            for j,((lo,hi),_) in enumerate(children):
                encloses(lo,add(p['lo'],mul(width,(F(j,16),F(0)))))
                encloses(hi,add(p['lo'],mul(width,(F(j+1,16),F(0)))))
        elif p['synthetic']:
            encloses(children[0][0][0],p['lo']);encloses(children[-1][0][1],p['hi'])
        else:assert children[0][0][0]==p['lo'] and children[-1][0][1]==p['hi']
        for ((_,hi),_),((lo,_),_) in zip(children,children[1:]):assert hi==lo
        free_children+=len(children)
        for (_, _),etas in children:
            etas.sort(key=lambda r:r[0][0]);encloses(etas[0][0],(F(0),F(0)));encloses(etas[-1][1],(F(1),F(0)))
            for (_,hi),(lo,_) in zip(etas,etas[1:]):assert hi==lo
            if method[p['chart']]==1 and p['dep']==0:
                assert all(lo[1]==hi[1]==0 and 0<hi[0]-lo[0]<=F(1,4) for lo,hi in etas)
            else:
                n=len(etas)
                for j,(lo,hi) in enumerate(etas):encloses(lo,(F(j,n),F(0)));encloses(hi,(F(j+1,n),F(0)))
    for slot,r in ledger.items():
        assert r['plan_complete']=='1' and r['failed_chart']=='-1'
        assert counts[slot]==int(r['new_regular_panels']) and old_oblique[slot]==int(r['reused_oblique_panels'])
        n=counts[slot]+old_oblique[slot];assert n==int(r['combined_panels']) and n*256==int(r['planned_positive_queries']) and n<=51072
    assert c['completed_typed_plans']==29 and c['new_complex_axis_field_queries']==c['new_event_root_jobs']==c['new_inverse_solves']==c['new_physical_quadrature_queries']==c['geometry_reconstructions']==0
    previous=json.loads((ROOT/'shared_atlas_adoption_v5/contract.json').read_text());full=[]
    for r in previous['target_cost']:
        target=int(r['target']);pending=[x for x in ledger.values() if int(x['target'])==target]
        assert len(pending)==r['incomplete_atlas_slots']
        new=r['active_new_panels']+sum(int(x['combined_panels']) for x in pending);retained=r['retained_original_target_panels'];queries=(new+retained)*256
        assert queries<=13074432
        full.append(dict(target=target,complete_candidate_plan_panels=new+retained,planned_positive_queries=queries,target_query_limit=13074432,cost_pass=True,geometry_density_clearance_status="not_decided_by_cost_audit",all12_analytic_remainders_complete=False))
    assert len(full)==24 and sum(r['complete_candidate_plan_panels'] for r in full)==34955+1047+50782+4384==91168
    print(json.dumps(dict(scope='independent exact full205-parent/50782-panel coverage and124 convex-chart/origin metrics; retained4384 oblique plan; complete24 candidate costs; not complex clearance or all12 physics',
        run=str(run),independent_pending29_typed_plan_audit_pass=True,
        candidate_plans_complete=29,new_regular_panels=50782,reused_oblique_panels=4384,
        actual_parent_rectangles=205,graded_free_children=free_children,chart_metric_records=124,
        panel_policies=dict(policy),all24_complete_candidate_costs_known=True,
        complete_candidate_panels=91168,complete_candidate_queries=91168*256,
        maximum_target_queries=max(r['planned_positive_queries'] for r in full),target_cost=full,
        remaining_incomplete_geometry_density_atlas_slots=29,all12_analytic_remainders_complete=False,full_quadrature_allowed=False),indent=2))

if __name__=='__main__':main()

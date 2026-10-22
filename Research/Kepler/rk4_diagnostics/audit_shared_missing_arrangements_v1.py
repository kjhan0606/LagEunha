"""Exact retained canonical root transforms and recovered88 geometry accounting.

Internal scalar-root and real graph enclosure premises remain those of the
original native implementation. This is not an analytic quadrature certificate.
"""
from collections import defaultdict,Counter
from fractions import Fraction as F
from pathlib import Path
import json,struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball,add,mul,neg,div,encloses,verify,string
from audit_shared_regular_q_transfer_v1 import rows,Q,sha

ROOT=Path(__file__).resolve().parent

def main():
    run=ROOT/'shared_missing_arrangements_v1_mpi1';c=verify(run)
    inventory={(int(r['slot_id']),int(r['chart'])):r for r in rows(ROOT/'shared_remaining29_recovery_inventory_v1/actual_chart_recovery_inventory.csv') if r['retained_old_kind']=='completed_clipped_arrangement_metadata_missing'}
    assert len(inventory)==88
    hc,charts=wire(run/'recovered_chart_snapshots_rank0.bin');hs,strips=wire(run/'recovered_actual_strips_rank0.bin');hf,families=wire(run/'boundary_family_bindings_rank0.bin')
    assert (hc[3],hs[3],hf[3])==(1520,128,288) and len(charts)==88 and len(strips)==169
    state_path=ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin';assert sha(state_path)=='74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005'
    state=state_path.read_bytes();fm={};sb=defaultdict(list);fb=defaultdict(list)
    for z in charts:
        slot,source,target,chart,complete,empty,whole,mask,nstrip,windows=struct.unpack_from('<10Q',z)
        key=slot,chart;assert key in inventory and key not in fm;fm[key]=z;r=inventory[key]
        assert (source,target)==(int(r['source']),int(r['target'])) and complete==1 and empty==whole==0
        dep,trim=struct.unpack_from('<2i',z,80);assert dep==int(r['old_dependent']) and nstrip==int(r['old_strips'])
        apex,part=struct.unpack_from('<2i',z,96);assert chart==2*apex+part
        cut=quad(z,432);assert cut==1 and trim==-1 # Actual recovered result, not an assumed default.
        A=[F(int(j==apex)) for j in range(3)]
        B=[F(1,3)]*3 if part else [F(int(j in (apex,(apex+1)%3)),2) for j in range(3)]
        C=[F(int(j in (apex,(apex+2)%3)),2) for j in range(3)] if part else [F(1,3)]*3
        for a,node in enumerate((A,B,C)):
            for j in range(3):encloses(ball(z,112+32*(3*a+j)),(node[j],F(0)))
        assert ball(z,400)[0]>ball(z,400)[1]>=0
        for j in range(23):
            if mask&(1<<j):assert abs(ball(z,592+32*j)[0])>ball(z,592+32*j)[1]
    assert set(fm)==set(inventory)
    canonical={}
    for k in range(4):
        for r in rows(ROOT/f'shared_registry_roots_v1_mpi4/canonical_event_root_packets_rank{k}.csv'):
            key=int(r['canonical_edge']),int(r['first_target']),int(r['second_target']);assert r['complete']=='1'
            canonical.setdefault(key,[])
            if r['root_index']!='-1':canonical[key].append((Q(r['t']),Q(r['t_radius'])))
    kinds=Counter();seen_family=set();internal={};family_calls=[]
    for b in families:
        slot,source,target,chart,face,side,kind,complete,windows,index,count=struct.unpack_from('<11Q',b)
        key=slot,chart;z=fm[key];assert (source,target)==struct.unpack_from('<2Q',z,8)
        assert complete and 0<=face<23 and side<2 and struct.unpack_from('<Q',z,56)[0]&(1<<face)
        call=slot,chart,face,side;fb[call].append(b)
        if index in (0,(1<<64)-1):
            assert call not in seen_family;seen_family.add(call);kinds[kind]+=1;family_calls.append(call)
        if count:assert 0<=index<count and 0<ball(b,96)[0]-ball(b,96)[1] and ball(b,96)[0]+ball(b,96)[1]<1
        else:assert index==(1<<64)-1
    for call,records in fb.items():
        z=fm[call[:2]];count=struct.unpack_from('<Q',records[0],80)[0]
        assert len(records)==max(1,count)
        assert [struct.unpack_from('<Q',b,72)[0] for b in records]==(list(range(count)) if count else [(1<<64)-1])
        kind=struct.unpack_from('<Q',records[0],48)[0];identity=string(records[0],128,160)
        if kind==1:
            tag,edge,first,second,cuttext=identity.split(';');assert tag=='canonical'
            cut=Q(cuttext.split('=',1)[1]);assert cut==quad(z,432)
            source=struct.unpack_from('<Q',z,8)[0];at=64+1064*source;apex,part=struct.unpack_from('<2i',z,96)
            side=call[3];dep=struct.unpack_from('<i',z,80)[0];assert dep==1 and side==int(bool(part))
            edge=int(edge);local=next(j for j in range(3) if struct.unpack_from('<i',state,at+344+200*j+184)[0]==edge)
            start,end=struct.unpack_from('<2i',state,at+944+8*local);other=(apex+(2 if part else 1))%3
            assert {start,end}=={apex,other}
            expected=[]
            for t in canonical[edge,int(first),int(second)]:
                u=div(mul((F(2),F(0)),t if start==apex else add((F(1),F(0)),neg(t))),(cut,F(0)))
                if u[0]+u[1]<0 or u[0]-u[1]>1:continue
                assert u[0]-u[1]>0 and u[0]+u[1]<1;expected.append(u)
            expected.sort(key=lambda x:x[0]);assert len(expected)==count
            for b,u in zip(records,expected):encloses(ball(b,96),u)
        elif kind in (3,4,5):
            apex,part=struct.unpack_from('<2i',z,96);dep=struct.unpack_from('<i',z,80)[0]
            _,source,target,competitor,label,cuttext=identity.split(';')
            geometry=apex if dep==1 else 3+apex if call[3]==0 else 6+call[1]
            assert label=='line='+str(geometry) and Q(cuttext.split('=',1)[1])==quad(z,432)
            actual=[b[96:128] for b in records] if count else []
            if kind==4:assert identity in internal and internal[identity]==actual
            else:
                assert identity not in internal;internal[identity]=actual
                if kind==5:assert dep==0 and call[3]==0 and count==0 and struct.unpack_from('<Q',records[0],64)[0]==0
                else:assert 0<struct.unpack_from('<Q',records[0],64)[0]<=256
        else:assert kind==2 and count==0
    expected_calls={(slot,chart,j,b) for (slot,chart),z in fm.items() for j in range(23) if struct.unpack_from('<Q',z,56)[0]&(1<<j) for b in (0,1)}
    assert set(fb)==expected_calls
    for b in strips:
        slot,source,target,chart,index=struct.unpack_from('<5Q',b);key=slot,chart;z=fm[key]
        ci,dep,lo,hi=struct.unpack_from('<4i',b,48)
        assert (ci,dep)==(chart,struct.unpack_from('<i',z,80)[0]) and (source,target)==struct.unpack_from('<2Q',z,8)
        assert index==len(sb[key]);sb[key].append(b)
        assert lo==-1 or struct.unpack_from('<Q',z,56)[0]&(1<<lo)
        assert hi==-1 or struct.unpack_from('<Q',z,56)[0]&(1<<hi)
        assert ball(b,64)[0]+ball(b,64)[1]<ball(b,96)[0]-ball(b,96)[1]
        assert ball(b,64)[0]-ball(b,64)[1]>=0 and ball(b,96)[0]+ball(b,96)[1]<=1
    for key,z in fm.items():assert len(sb[key])==struct.unpack_from('<Q',z,64)[0]
    assert kinds=={1:149,3:124,4:52,5:27},kinds
    assert sum(struct.unpack_from('<Q',fb[call][0],64)[0] for call in family_calls if struct.unpack_from('<Q',fb[call][0],48)[0]==3)==c['new_internal_root_windows']==337
    assert c['completed_recovered_charts']==88 and c['stored_actual_strips']==169
    assert c['new_canonical_root_jobs']==c['new_inverse_solves']==c['new_physical_quadrature_queries']==0 and not c['full_quadrature_allowed']
    print(json.dumps(dict(scope='exact88 recovered frame/status/169strip bindings and complete active-side families; canonical parameter transforms and internal cache identities; retained native scalar-root/graph enclosure premises',
        independent_missing88_arrangement_audit_pass=True,run=str(run),recovered_charts=88,recovered_strips=169,
        reused_canonical_root_families=149,new_internal_root_families=124,new_internal_root_windows=337,
        reused_internal_root_families=52,constant_apex_no_root_proofs=27,
        whole_or_excluded72_chart_reference_snapshots_still_unrecorded=True,
        remaining_incomplete_atlas_slots=29,all12_analytic_remainders_complete=False,full_quadrature_allowed=False),indent=2))

if __name__=='__main__':main()

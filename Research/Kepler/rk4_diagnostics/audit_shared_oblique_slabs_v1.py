"""Exact saved actual-event coverage and selected-branch slab accounting."""
from collections import defaultdict,Counter
from fractions import Fraction as F
from pathlib import Path
import json
import struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball,add,neg,mul,div,encloses,verify,string
from audit_shared_regular_q_transfer_v1 import rows,sha
ROOT=Path(__file__).resolve().parent
def sub(a,b):return add(a,neg(b))
def scale(a,c):return mul(a,(c,0))
def positive(a):return a[0]>a[1]
def negative(a):return a[0]<-a[1]
def coordinates(v,apex,part,X,Y):
    ax,ay=v[0];bx,by=v[1];cx,cy=v[2];det=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax)
    x=sub(X,(ax,0));y=sub(Y,(ay,0));l1=add(scale(x,(cy-ay)/det),scale(y,-(cx-ax)/det));l2=add(scale(x,-(by-ay)/det),scale(y,(bx-ax)/det));l=[sub((1,0),add(l1,l2)),l1,l2]
    B=[F(1,3)]*3 if part else [F(int(k in (apex,(apex+1)%3)),2) for k in range(3)]
    C=[F(int(k in (apex,(apex+2)%3)),2) for k in range(3)] if part else [F(1,3)]*3
    i,j=(apex+1)%3,(apex+2)%3;det0=B[i]*C[j]-B[j]*C[i];assert det0>0
    u=scale(sub(scale(l[i],C[j]),scale(l[j],C[i])),1/det0);q=scale(sub(scale(l[j],B[i]),scale(l[i],B[j])),1/det0);s=add(u,q)
    return s,div(q,s)
def main():
    run=ROOT/'shared_oblique_slabs_v1_mpi1';c=verify(run);previous=ROOT/'shared_pending29_rebind_v1_mpi1';inside=ROOT/'shared_oblique_internal_events_v1_mpi1'
    # The original runner verified the entire parent output manifest before
    # launch; it includes the consumed node packet even though that path was
    # omitted from the redundant direct-input list. Verify the same chain.
    for line in (previous/'output.sha256').read_text().splitlines():
        digest,path=line.split('  ',1);assert sha(Path(path))==digest
    _,frames0=wire(previous/'oblique_coordinate_frames_rank0.bin');_,frames1=wire(run/'unchanged_oblique_frames_rank0.bin');assert frames1==frames0
    frames={struct.unpack_from('<Q',b)[0]:b for b in frames0};assert len(frames)==14
    _,nodes0=wire(previous/'actual_reference_events_rank0.bin');nodes={(struct.unpack_from('<Q',b)[0],struct.unpack_from('<Q',b,24)[0]):b for b in nodes0}
    _,roots0=wire(inside/'internal_chart_boundary_events_rank0.bin');roots={}
    for b in roots0:
        slot,source,target,chart,face,side,index,count,complete,windows=struct.unpack_from('<10Q',b);assert complete
        if index!=(1<<64)-1:roots[slot,side,face,index]=b
    state=(ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin').read_bytes();assert sha(ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin')=='74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005'
    vertices={s:[tuple(F(x) for x in struct.unpack_from('<2d',state,64+s*1064+16*k)) for k in range(3)] for s in {struct.unpack_from('<Q',b,8)[0] for b in frames0}}
    he,events=wire(run/'ordered_projection_events_rank0.bin');hs,slabs=wire(run/'actual_oblique_slabs_rank0.bin');hp,proofs=wire(run/'slab_face_proofs_rank0.bin');assert (he[3],hs[3],hp[3])==(400,304,176)
    eb=defaultdict(list);sb=defaultdict(list);pb=defaultdict(list);kinds=Counter();seenroots=set();seennodes=set();corners=set()
    for b in events:
        slot,source,target,chart,index,kind=struct.unpack_from('<6Q',b);face,side,node=struct.unpack_from('<3q',b,48);f=frames[slot];assert (source,target,chart)==struct.unpack_from('<3Q',f,8);assert index==len(eb[slot]);eb[slot].append(b);kinds[kind]+=1;identity=string(b,208,192)
        if kind==0:
            key=slot,side;assert key not in corners;corners.add(key);assert side in range(4) and face==node==-1 and b[80:144]==f[864+64*side:864+64*(side+1)]
            assert ball(b,144)==(quad(f,1856) if side in (0,3) else F(1),F(0)) and ball(b,176)==(F(side>=2),F(0))
        elif kind==1:
            label=identity.split(';');assert label[:2]==['internal',str(slot)];key=(slot,side,face,int(label[4]));assert key in roots and key not in seenroots;seenroots.add(key);r=roots[key];assert b[80:144]==r[240:304] and b[144:208]==r[112:176]
        else:
            key=slot,node;assert key in nodes and key not in seennodes;seennodes.add(key);r=nodes[key];assert identity==string(r,288,160) and struct.unpack_from('<Q',r,40)[0]
            assert kind==(3 if struct.unpack_from('<Q',r,32)[0]==2 else 2)
    assert corners=={(s,k) for s in frames for k in range(4)} and seenroots==set(roots)
    skipped=Counter();included=Counter()
    for key,r in nodes.items():
        slot,node=key
        if slot not in frames:continue
        f=frames[slot];source=struct.unpack_from('<Q',f,8)[0];apex,part=struct.unpack_from('<2i',f,48);X,Y=ball(r,128),ball(r,160);cut=quad(f,1856)
        if all(ball(r,128+32*k)==ball(f,384+32*k) and ball(r,128+32*k)[1]==0 for k in range(2)):
            assert key not in seennodes;skipped[slot]+=1;continue
        s,t=coordinates(vertices[source],apex,part,X,Y);other=(apex+(2 if part else 1))%3;a=64+source*1064;edgekey=None
        for j in range(3):
            x,y=struct.unpack_from('<2i',state,a+944+8*j)
            if {x,y}=={apex,other}:edgekey=struct.unpack_from('<i',state,a+344+200*j+184)[0]
        assert edgekey is not None;canonical=struct.unpack_from('<Q',r,32)[0]==1 and string(r,288,160).startswith(f'e;{edgekey};')
        if canonical:t=(F(bool(part)),F(0))
        outside=s[0]+s[1]<cut or s[0]-s[1]>1 or t[0]+t[1]<0 or t[0]-t[1]>1
        if key in seennodes:
            assert not outside and s[0]-s[1]>cut and s[0]+s[1]<1 and (canonical or 0<t[0]-t[1]<=t[0]+t[1]<1);included[slot]+=1
            b=next(b for b in eb[slot] if struct.unpack_from('<q',b,64)[0]==node);encloses(ball(b,144),s);encloses(ball(b,176),t)
        else:assert outside;skipped[slot]+=1
    for b in slabs:
        slot,source,target,chart,index,empty=struct.unpack_from('<6Q',b);assert index==len(sb[slot]);sb[slot].append(b);f=frames[slot];geolo,geohi,lower,upper=struct.unpack_from('<4q',b,48);assert geolo in range(4) and geohi in range(4)
        assert b[80:112]==eb[slot][index][80:112] and b[112:144]==eb[slot][index+1][80:112];sample=ball(b,144);assert sample[1]==0 and ball(b,80)[0]+ball(b,80)[1]<sample[0]<ball(b,112)[0]-ball(b,112)[1]
        hit=[]
        for side in range(4):
            nxt=(side+1)%4;u0,v0=ball(f,864+64*side),ball(f,896+64*side);u1,v1=ball(f,864+64*nxt),ball(f,896+64*nxt)
            if min(u0[0],u1[0])<sample[0]<max(u0[0],u1[0]):hit.append((add(v0,div(mul(sub(sample,u0),sub(v1,v0)),sub(u1,u0))),side))
        assert len(hit)==2;hit.sort(key=lambda x:x[0][0]);assert [h[1] for h in hit]==[geolo,geohi];encloses(ball(b,176),hit[0][0]);encloses(ball(b,208),hit[1][0]);assert ball(b,176)[0]+ball(b,176)[1]<ball(b,208)[0]-ball(b,208)[1]
    for b in proofs:
        slot,slab,face,case,iterations=struct.unpack_from('<5Q',b);f=frames[slot];assert struct.unpack_from('<Q',f,40)[0]&(1<<face);assert b[48:80]==f[1120+32*face:1152+32*face];left,right,root=ball(b,80),ball(b,112),ball(b,144)
        if case==0:assert negative(left) and negative(right) and iterations==0
        elif case==1:assert positive(left) and positive(right) and iterations==0
        else:assert case==2 and ((negative(left) and positive(right)) or (positive(left) and negative(right))) and 0<iterations<=256 and root[1]>=0
        pb[slot,slab].append(b)
    empty_count=0
    for slot,sl in sb.items():
        mask=struct.unpack_from('<Q',frames[slot],40)[0];active=[i for i in range(23) if mask&(1<<i)]
        for j,b in enumerate(sl):
            ps=pb[slot,j];assert [struct.unpack_from('<Q',p,16)[0] for p in ps]==active
            lo,hi=ball(b,176),ball(b,208);lower=upper=-1;empty=False
            for p in ps:
                face=struct.unpack_from('<Q',p,16)[0];case=struct.unpack_from('<Q',p,24)[0]
                if case==1:empty=True
                if case!=2:continue
                root=ball(p,144);derivative=ball(p,48)
                if positive(derivative):
                    if root[0]+root[1]<hi[0]-hi[1]:hi=root;upper=face
                    else:assert root[0]-root[1]>hi[0]+hi[1]
                else:
                    assert negative(derivative)
                    if root[0]-root[1]>lo[0]+lo[1]:lo=root;lower=face
                    else:assert root[0]+root[1]<lo[0]-lo[1]
            assert ball(b,240)==lo and ball(b,272)==hi and struct.unpack_from('<2q',b,64)==(lower,upper)
            if lo[0]-lo[1]>hi[0]+hi[1]:empty=True
            if not empty:assert lo[0]+lo[1]<hi[0]-hi[1]
            assert bool(struct.unpack_from('<Q',b,40)[0])==empty;empty_count+=empty
    for r in rows(run/'actual_oblique_arrangements.csv'):
        slot=int(r['slot_id']);ev=eb[slot];assert r['ordered']==r['complete']=='1' and len(ev)==int(r['events']) and len(sb[slot])==len(ev)-1==int(r['slabs']) and skipped[slot]==int(r['skipped_outside_original_chart_events'])
        for a,b in zip(ev,ev[1:]):assert ball(a,80)[0]+ball(a,80)[1]<ball(b,80)[0]-ball(b,80)[1]
    assert len(events)==c['saved_projection_events']==112 and len(slabs)==c['saved_slabs']==98 and len(proofs)==c['saved_face_sample_proofs']==294 and empty_count==c['explicit_empty_slabs']==58 and c['completed_oblique_chart_arrangements']==14
    assert all(c[k]==0 for k in ('new_inverse_solves','new_canonical_root_jobs','new_internal_boundary_root_jobs','normal_field_or_panel_replays','new_physical_quadrature_queries')) and not c['full_quadrature_allowed']
    print(json.dumps(dict(scope='exact original corner/internal/canonical/power-event inclusion/exclusion, strict ordering/shared endpoint coverage and selected-branch/empty-slab proof accounting; retains native field/root and original topology premises',
        run=str(run),independent_actual_oblique_slab_accounting_pass=True,completed_oblique_chart_arrangements=14,projection_events=112,slabs=98,nonempty_slabs=40,empty_slabs=58,face_sample_proofs=294,
        event_kind_counts=dict(kinds),typed_panels_and_complex_density_clearance_still_pending=True,remaining_incomplete_atlas_slots=29,full_quadrature_allowed=False,
        inherited_reference_event_input_checked_through_full_parent_output_manifest=True),indent=2))
if __name__=='__main__':main()

"""Frozen new regular-axis prefix/full-field/boundary/branch accounting.

Retains the original native analytic-field/libm/implicit-domain premises.
No physical query or time evolution is performed by this audit.
"""
from collections import defaultdict,Counter
from fractions import Fraction as F
from pathlib import Path
import json,struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball,add,mul,neg,encloses,verify
from audit_shared_regular_q_transfer_v1 import rows,sha

ROOT=Path(__file__).resolve().parent
MAX=(1<<64)-1

def minus(a,b):return add(a,neg(b))

def main(run_name='shared_pending_clearance_v1_mpi1',
         plan_name='shared_pending_atlas_v1_mpi1',
         panel_filename='new_regular_typed_panels_rank0.bin',
         expected_panels=50782, expected_slots=29, verification=verify):
    run=ROOT/run_name;c=verification(run)
    plan=ROOT/plan_name
    proof=json.loads((ROOT/'SHARED_PENDING_TYPED_ATLAS_20261021.json').read_text());assert proof['independent_pending29_typed_plan_audit_pass']
    approved=json.loads((ROOT/'shared_atlas_adoption_v5/contract.json').read_text())['required_proof_helpers']
    for original,digest in approved.items():assert sha(run/Path(original).name)==digest
    _,panels=wire(plan/panel_filename)
    pm={(struct.unpack_from('<Q',b)[0],struct.unpack_from('<Q',b,16)[0]):b for b in panels}
    counts=Counter(slot for slot,index in pm);assert len(pm)==expected_panels and len(counts)==expected_slots
    ha,axes=wire(run/'regular_axis_proofs_rank0.bin');hp,points=wire(run/'passed_regular_full_fields_rank0.bin');hb,boundaries=wire(run/'prepared_regular_boundaries_rank0.bin')
    assert (ha[3],hp[3],hb[3])==(624,1680,1552)
    byslot=defaultdict(list);bb=defaultdict(list)
    for b in boundaries:
        slot,source,target,ch,index,axis=struct.unpack_from('<6Q',b);p=pm[slot,index]
        assert (source,target)==(struct.unpack_from('<Q',p,8)[0],struct.unpack_from('<Q',p,24)[0])
        assert ch==struct.unpack_from('<i',p,48)[0] and axis<2 and not p[192]
        bb[slot].append(b)
        dep=struct.unpack_from('<i',p,36)[0]
        lower,upper=struct.unpack_from('<2i',p,40)
        offsets=(560,1168) if dep==1 else (1296,1424)
        for face,at in zip((lower,upper),offsets):
            if face<0:continue
            q,inv,R,prior=quad(b,at+32),quad(b,at+48),quad(b,at+80),quad(b,at+96)
            assert 0<=q<1 and 0<=inv<1 and 0<R<=prior
        if dep==1:
            assert b[80+224]==b[80+464]==b[688+224]==b[688+464]==0
    point_index=0;quadratic=0;state=(ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin').read_bytes()
    for a in axes:
        slot,source,target,index,axis,passed,boundary=struct.unpack_from('<7Q',a);p=pm[slot,index]
        assert (source,target)==(struct.unpack_from('<Q',p,8)[0],struct.unpack_from('<Q',p,24)[0]) and axis<2
        byslot[slot].append(a)
        U=(F(1,2),F(1,2) if axis else F(5,4));E=(F(1,2),F(5,4) if axis else F(1,2))
        if p[192]:
            sl,sh,tl,th=ball(p,240),ball(p,272),ball(p,304),ball(p,336)
            encloses(ball(a,64),add(sl,mul(minus(sh,sl),U)));encloses(ball(a,96),add(tl,mul(minus(th,tl),E)))
        else:
            encloses(ball(a,64),add(ball(p,64),mul(minus(ball(p,96),ball(p,64)),U)))
            encloses(ball(a,96),add(ball(p,128),mul(minus(ball(p,160),ball(p,128)),E)))
        if not passed:continue
        f=points[point_index];point_index+=1
        assert a[:40]==f[:40] and a[224:256]==f[1440:1472] and a[192:224]==f[1472:1504] and a[256:288]==f[432:464]
        assert ball(a,256)[0]>ball(a,256)[1]>=0
        assert 0<=quad(a,320)<1 and 0<=quad(a,336)<1 and 0<=quad(a,368)<=quad(a,384)
        if p[192]:
            assert boundary==MAX and a[64:96]==a[128:160] and a[96:128]==a[160:192]
            expected=mul(mul(ball(f,1360),minus(ball(p,272),ball(p,240))),minus(ball(p,336),ball(p,304)))
        else:
            dep=struct.unpack_from('<i',p,36)[0];lo,hi=struct.unpack_from('<2i',p,40)
            free_width=minus(ball(p,96),ball(p,64));eta_width=minus(ball(p,160),ball(p,128))
            if lo>=0 or hi>=0:
                assert boundary<len(bb[slot]);b=bb[slot][boundary]
                assert struct.unpack_from('<Q',b,40)[0]==axis and b[48:80]==a[64:96]
                first=struct.unpack_from('<Q',b,32)[0];assert first<=index and pm[slot,first][48:128]==p[48:128]
                if dep==1:
                    width=add(add(mul(ball(a,64),minus(ball(b,688),ball(b,80))),minus(ball(b,720),ball(b,112))),minus(ball(b,752),ball(b,144)))
                    expected=mul(mul(mul(ball(p,672),width),free_width),eta_width)
                else:
                    width=minus(ball(b,1424),ball(b,1296))
                    expected=mul(mul(mul(ball(f,1360),width),free_width),eta_width)
            else:
                assert boundary==MAX
                expected=mul(mul(ball(f,1360),free_width),eta_width)
        encloses(ball(a,192),expected)
        method=struct.unpack_from('<Q',f,1504)[0]
        if method:
            assert method==1;quadratic+=1;apex,part=struct.unpack_from('<2i',p,368)
            v=[tuple(F(x) for x in struct.unpack_from('<2d',state,64+1064*source+16*j)) for j in range(3)]
            assert v[apex]==(0,0)
            B=[F(1,3)]*3 if part else [F(int(j in (apex,(apex+1)%3)),2) for j in range(3)]
            C=[F(int(j in (apex,(apex+2)%3)),2) for j in range(3)] if part else [F(1,3)]*3
            u=[sum(B[j]*v[j][k] for j in range(3)) for k in range(2)]
            d=[sum((C[j]-B[j])*v[j][k] for j in range(3)) for k in range(2)]
            aa=sum(x*x for x in u);ab=2*sum(u[k]*d[k] for k in range(2));ac=sum(x*x for x in d)
            for j,value in enumerate((aa,ab,ac)):encloses(ball(f,1520+32*j),(value,F(0)))
            t=ball(a,160);distance,gap,sc,derivative=(quad(f,1616+16*j) for j in range(4))
            alpha=-ab/(2*ac);beta2=(4*aa*ac-ab*ab)/(4*ac*ac)
            assert beta2>0 and 0<gap<=distance-t[1] and distance*distance<=(t[0]-alpha)**2+beta2
            assert sc>0 and sc*sc<=ac and derivative>=abs(ab+2*ac*t[0])+2*ac*t[1]
    assert point_index==len(points)==c['saved_full_fields']==c['new_complex_axis_passes']
    ledger={int(r['slot_id']):r for r in rows(run/'regular_geometry_density_clearance.csv')};assert set(ledger)==set(counts)
    complete=0;first=[]
    for slot,n in counts.items():
        r=ledger[slot];ar=byslot[slot];expected=[(i,j) for i in range(n) for j in (0,1)]
        assert [(struct.unpack_from('<Q',a,24)[0],struct.unpack_from('<Q',a,32)[0]) for a in ar]==expected[:len(ar)]
        passes=sum(struct.unpack_from('<Q',a,40)[0] for a in ar)
        assert n==int(r['regular_panels']) and len(ar)==int(r['checked_axes']) and passes==int(r['passed_axes'])
        assert len(bb[slot])==int(r['cached_boundary_records'])
        if r['clearance_complete']=='1':assert passes==len(ar)==2*n;complete+=1
        else:
            if ar:
                assert len(ar)==passes+1 and not struct.unpack_from('<Q',ar[-1],40)[0]
                first.append(dict(slot=slot,panel=struct.unpack_from('<Q',ar[-1],24)[0],axis=struct.unpack_from('<Q',ar[-1],32)[0],reason=r['reason']))
            else:
                assert passes==0 and r['reason']
                first.append(dict(slot=slot,panel=None,axis=None,reason=r['reason'],setup_failed_before_axis_queries=True))
    assert complete==c['completed_regular_slot_clearances'] and len(axes)==c['new_complex_axis_evaluations']
    assert c['previous_successful_field_replays']==c['new_event_root_jobs']==c['new_inverse_solves']==c['new_physical_quadrature_queries']==0 and c['state_unchanged'] and not c['full_quadrature_allowed']
    result=dict(scope='exact new regular-axis prefix/input/full-field/J/radius/measure and boundary-cache bindings; original positive-origin sqrt branch certificates; retains native analytic/implicit/libm premises, not all12 or trajectory',
        run=str(run),independent_pending_regular_clearance_audit_pass=True,
        new_axis_checks=len(axes),new_axis_passes=len(points),full_fields_saved=len(points),
        completed_regular_slot_clearances=complete,first_failures=first,
        prepared_boundary_records=len(boundaries),origin_quadratic_branch_certificates=quadratic,
        scheduled_regular_slot_clearances=expected_slots,
        all_scheduled_regular_clearance_complete=complete==expected_slots,
        all12_analytic_remainders_complete=False,full_quadrature_allowed=False,long_time_stability_demonstrated=False)
    if expected_slots==29:result['all29_new_regular_clearance_complete']=complete==29
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()

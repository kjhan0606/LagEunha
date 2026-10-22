"""Frozen limited internal-root accounting and exact coordinate/projection checks."""
from collections import defaultdict
from pathlib import Path
import json
import struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball,add,neg,mul,div,encloses,verify
from audit_shared_regular_q_transfer_v1 import rows
ROOT=Path(__file__).resolve().parent
def main():
    run=ROOT/'shared_oblique_internal_events_v1_mpi1';c=verify(run)
    _,original=wire(ROOT/'shared_pending29_rebind_v1_mpi1/oblique_coordinate_frames_rank0.bin');_,saved=wire(run/'unchanged_oblique_frames_rank0.bin');assert saved==original and len(saved)==14
    frames={struct.unpack_from('<Q',b)[0]:b for b in saved};h,events=wire(run/'internal_chart_boundary_events_rank0.bin');assert h[3]==496
    groups=defaultdict(list);rootcount=windows=0
    for b in events:
        slot,source,target,chart,face,side,index,count,complete,work=struct.unpack_from('<10Q',b);f=frames[slot];fs,source0,target0,chart0,passed,mask=struct.unpack_from('<6Q',f)
        assert (source,target,chart)==(source0,target0,chart0) and mask&(1<<face) and passed
        part=struct.unpack_from('<i',f,52)[0];assert side in ((1,2) if part==0 else (0,1)) # canonical side never replayed
        groups[slot,side,face].append(b)
        if index==(1<<64)-1:assert count==0;continue
        assert index<count;rootcount+=1;parameter=ball(b,80);assert 0<=parameter[0]-parameter[1]<=parameter[0]+parameter[1]<=1
        cut=(quad(f,1856),0)
        if side==0:s=add(cut,mul(add((1,0),neg(cut)),parameter));t=(0,0)
        elif side==1:s=(1,0);t=parameter
        else:s=add((1,0),neg(mul(add((1,0),neg(cut)),parameter)));t=(1,0)
        encloses(ball(b,112),s);encloses(ball(b,144),t)
        for k in range(2):
            A=ball(f,384+32*k);B=ball(f,608+64+32*k);C=ball(f,608+128+32*k)
            X=add(A,mul(ball(b,112),add(add(B,neg(A)),mul(ball(b,144),add(C,neg(B))))));encloses(ball(b,176+32*k),X)
        delta=[add(ball(b,176+32*k),neg(ball(f,384+32*k))) for k in range(2)]
        for axis,at in ((0,512),(1,448)):
            v=div(add(mul(ball(f,at),delta[0]),mul(ball(f,at+32),delta[1])),ball(f,576));encloses(ball(b,240+32*axis),v)
    completed=0;observed_jobs=0
    for r in rows(run/'internal_chart_boundary_status.csv'):
        slot=int(r['slot_id']);frame=frames[slot];mask=struct.unpack_from('<Q',frame,40)[0];part=struct.unpack_from('<i',frame,52)[0]
        expected=[(slot,s,face) for s in ((1,2) if part==0 else (0,1)) for face in range(23) if mask&(1<<face)]
        got=[k for k in groups if k[0]==slot];assert got==expected[:len(got)] and len(got)==int(r['new_internal_root_jobs']);observed_jobs+=len(got)
        for key in got:
            records=groups[key];metadata=[struct.unpack_from('<10Q',b) for b in records];first=metadata[0];count=first[7];work=first[9];windows+=work
            assert len(records)==max(1,count) and all(m[7:]==first[7:] for m in metadata)
            if count:assert [m[6] for m in metadata]==list(range(count))
            assert work<=257 # original guard counts the rejected257th window, if any
        if r['complete']=='1':assert got==expected and all(struct.unpack_from('<Q',b,64)[0] for k in got for b in groups[k]);completed+=1
    assert observed_jobs==c['new_internal_boundary_root_jobs'] and windows==c['new_internal_root_work_windows'] and rootcount==c['stored_internal_roots'] and len(events)==c['stored_rootset_rows'] and completed==c['completed_internal_boundary_root_charts']
    assert c['new_canonical_root_jobs']==c['new_inverse_solves']==c['normal_field_or_panel_replays']==c['new_physical_quadrature_queries']==0 and not c['full_quadrature_allowed'] and c['state_unchanged']
    print(json.dumps(dict(scope='same frozen native root guard retained; independent original14/two-internal-side ownership, complete/prefix rootset accounting and exact root coordinate/projection enclosures',
        run=str(run),internal_event_packet_accounting_pass=True,completed_internal_boundary_charts=completed,
        pending_internal_boundary_charts=14-completed,new_internal_root_jobs=observed_jobs,stored_internal_roots=rootcount,
        full_actual_event_slab_atlas_still_pending=True,remaining_incomplete_atlas_slots=29,full_quadrature_allowed=False),indent=2))
if __name__=='__main__':main()

"""Assemble explicit real geometry packets for174 charts; no new evaluator call."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import csv,json,struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball,encloses,verify
from audit_shared_regular_q_transfer_v1 import rows,sha

ROOT=Path(__file__).resolve().parent

def main():
    run=ROOT/'shared_chart72_metadata_v1_mpi1';c=verify(run)
    p88=ROOT/'shared_missing_arrangements_v1_mpi1'
    evidence=ROOT/'SHARED_MISSING_ARRANGEMENTS_20261021.json'
    a=json.loads(evidence.read_text());assert a['independent_missing88_arrangement_audit_pass']
    inventory={(int(r['slot_id']),int(r['chart'])):r for r in rows(ROOT/'shared_remaining29_recovery_inventory_v1/actual_chart_recovery_inventory.csv')}
    requested={key:r for key,r in inventory.items() if r['retained_old_kind'] in ('whole','excluded')};assert len(requested)==72
    h,closed=wire(run/'whole_excluded_chart_snapshots_rank0.bin');assert h[3]==1520 and len(closed)==72
    bindings={};kinds=Counter()
    def put(key,kind,file,index,width):
        assert key in inventory and key not in bindings
        bindings[key]=(kind,str(file),index,width)
        kinds[kind]+=1
    for i,b in enumerate(closed):
        slot,source,target,chart,complete,empty,whole,mask,nstrip,windows=struct.unpack_from('<10Q',b)
        key=slot,chart;r=requested[key]
        assert complete and (source,target)==(int(r['source']),int(r['target']))
        assert (empty,whole)==(int(r['retained_old_kind']=='excluded'),int(r['retained_old_kind']=='whole')) and nstrip==windows==0
        assert struct.unpack_from('<i',b,80)[0]==-1
        apex,part=struct.unpack_from('<2i',b,96);assert chart==2*apex+part
        cut=quad(b,432);assert cut==1 and struct.unpack_from('<i',b,84)[0]==-1
        A=[F(int(j==apex)) for j in range(3)]
        B=[F(1,3)]*3 if part else [F(int(j in (apex,(apex+1)%3)),2) for j in range(3)]
        C=[F(int(j in (apex,(apex+2)%3)),2) for j in range(3)] if part else [F(1,3)]*3
        for n,node in enumerate((A,B,C)):
            for j in range(3):encloses(ball(b,112+32*(3*n+j)),(node[j],F(0)))
        assert ball(b,400)[0]>ball(b,400)[1]>=0
        put(key,r['retained_old_kind'],run/'whole_excluded_chart_snapshots_rank0.bin',i,1520)
    file88=p88/'recovered_chart_snapshots_rank0.bin';_,v88=wire(file88)
    for i,b in enumerate(v88):
        slot,source,target,chart,complete=struct.unpack_from('<5Q',b);assert complete
        key=slot,chart;assert inventory[key]['retained_old_kind']=='completed_clipped_arrangement_metadata_missing'
        put(key,'clipped',file88,i,1520)
    oblique=ROOT/'shared_oblique_collected_suffix_v1_mpi1/unchanged_frames_rank0.bin'
    _,frames=wire(oblique)
    for i,b in enumerate(frames):
        slot,source,target,chart,complete=struct.unpack_from('<5Q',b);assert complete
        key=slot,chart;assert inventory[key]['retained_old_kind']=='failed_oblique_replaced'
        put(key,'oblique',oblique,i,2064)
    assert set(bindings)==set(inventory) and len(bindings)==174 and kinds=={'whole':36,'excluded':36,'clipped':88,'oblique':14}
    slots={slot for slot,ci in bindings};assert len(slots)==29
    assert all({ci for s,ci in bindings if s==slot}==set(range(6)) for slot in slots)
    assert c['complete_reference_snapshots']==72 and c['actual_nontrivial_radial_trims']==0
    assert c['new_root_jobs']==c['new_inverse_solves']==c['new_physical_quadrature_queries']==0
    assert c['state_unchanged'] and not c['full_quadrature_allowed']
    dest=ROOT/'shared_pending29_chart_assembly_v1';dest.mkdir()
    with (dest/'actual174_chart_bindings.csv').open('w',newline='') as out:
        w=csv.writer(out);w.writerow(('slot_id','source','target','chart','kind','typed_packet','record','width'))
        for (slot,ci),binding in sorted(bindings.items()):
            r=inventory[slot,ci];w.writerow((slot,r['source'],r['target'],ci,*binding))
    result=dict(scope='explicit typed real reference geometry for all174 charts of pending29; audited88 arrangements, native72 reference/trim metadata, previously adopted14 oblique frames; not typed integration atlas or physics',
        independent72_reference_metadata_and174_binding_audit_pass=True,
        full174_real_geometry_packet_bindings_complete=True,slots=29,charts=174,chart_kinds=dict(kinds),
        new_root_jobs_for72=0,new_field_or_quadrature_queries_in_assembly=0,
        actual_nontrivial_radial_trims_in72=0,remaining_incomplete_atlas_slots=29,
        pending='origin15 and other regular chart typed panels/complex clearance, target cost, all12 remainder, evolved heat and actual native long evolution',
        all12_analytic_remainders_complete=False,full_quadrature_allowed=False,long_time_stability_demonstrated=False)
    (dest/'contract.json').write_text(json.dumps(result,indent=2)+'\n')
    inputs=[ROOT/'shared_remaining29_recovery_inventory_v1/actual_chart_recovery_inventory.csv',evidence,
        run/'whole_excluded_chart_snapshots_rank0.bin',file88,p88/'recovered_actual_strips_rank0.bin',oblique,
        ROOT/'SHARED_CHART72_METADATA_POLICY_20261021.md']
    (dest/'input.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in inputs))
    (dest/'source.sha256').write_text(sha(Path(__file__))+'  '+str(Path(__file__).resolve())+'\n')
    (dest/'output.sha256').write_text(''.join(sha(dest/n)+'  '+str(dest/n)+'\n' for n in ('contract.json','actual174_chart_bindings.csv')))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()

"""Remaining fixed runtime work; reuse saved implicit caches, not arrangements."""
from collections import defaultdict,Counter
from pathlib import Path
import csv,json,struct
from audit_shared_regular_q_transfer_v1 import rows,sha

ROOT=Path(__file__).resolve().parent
MAX=(1<<64)-1

def main():
    dest=ROOT/'shared_mixed_runtime_inputs_v1';dest.mkdir()
    registry=ROOT/'shared_atlas_adoption_v7';adopt=ROOT/'shared_native_transfer_adoption_v1/contract.json'
    a=json.loads(adopt.read_text());assert a['completed_native_source_overlap_records']==64
    remaining=set(a['remaining_target_ids']);packets={int(r['packet_id']):r for r in rows(registry/'typed_packet_registry.csv')}
    inputs=[adopt,registry/'contract.json',registry/'panel_evidence_registry.csv',registry/'typed_packet_registry.csv']
    raw={};headers={};boundary_slots={}
    def load(path):
        path=str(path)
        if path not in raw:
            data=Path(path).read_bytes();h=struct.unpack_from('<8Q',data)
            assert h[0] in (0x4753415749544e31,0x475341544c415331) and h[1:3]==(1,1) and h[6:]==(24,384) and len(data)==64+h[3]*h[4]
            raw[path]=data;headers[path]=h;inputs.append(Path(path))
        return raw[path],headers[path]
    def record(path,index):
        data,h=load(path);assert 0<=index<h[4]
        return memoryview(data)[64+h[3]*index:64+h[3]*(index+1)]
    def registered(pid,index):
        p=packets[pid];data,h=load(p['path']);assert h[3]==int(p['record_width']) and h[4]==int(p['records'])
        return record(p['path'],index)
    def boundary(path,slot,local):
        path=str(path);data,h=load(path);assert h[3]==1552
        if path not in boundary_slots:
            byslot=defaultdict(list)
            for j in range(h[4]):byslot[struct.unpack_from('<Q',data,64+1552*j)[0]].append(j)
            boundary_slots[path]=byslot
        return boundary_slots[path][slot][local]
    assembly=ROOT/'shared_pending29_chart_assembly_v1/actual174_chart_bindings.csv';inputs.append(assembly)
    frames={(int(r['slot_id']),int(r['chart'])):r for r in rows(assembly) if r['kind']=='oblique'}
    groups={};group_rows=[];jobs=[];classes=Counter();slots=set();old_naive_branches=0
    for r in rows(registry/'panel_evidence_registry.csv'):
        target=int(r['target'])
        if target not in remaining:continue
        slot,source,index=map(int,(r['slot_id'],r['source'],r['native_panel']));slots.add(slot)
        pid,rec=int(r['typed_packet_id']),int(r['typed_record']);p=registered(pid,rec)
        row=dict(slot_id=slot,source=source,target=target,panel_kind=r['panel_kind'],native_panel=index,
                 typed_packet_id=pid,typed_record=rec,geometry_group=-1,prepared_boundary_packet='',prepared_boundary_record=-1,
                 oblique_frame_packet='',oblique_frame_record=-1)
        if r['panel_kind']=='regular':
            assert len(p)==704 and struct.unpack_from('<4Q',p)==(slot,source,index,target)
            whole=bool(p[192]);dep,lo,hi=struct.unpack_from('<3i',p,36)
            classes['whole' if whole else 'angular' if dep==1 else 'radial']+=1
            if not whole and (lo>=0 or hi>=0):
                proof_pid,proof_index=map(int,(r['axis1_proof_packet'],r['axis1_proof_record']))
                proof=registered(proof_pid,proof_index);assert len(proof)==624
                s,src,tgt,idx,axis,passed,local=struct.unpack_from('<7Q',proof)
                assert (s,src,tgt,idx,axis,passed)==(slot,source,target,index,1,1) and local!=MAX
                path=Path(packets[proof_pid]['path']).parent/'prepared_regular_boundaries_rank0.bin'
                actual=boundary(path,slot,local);b=record(path,actual)
                assert struct.unpack_from('<4Q',b)==(slot,source,target,struct.unpack_from('<i',p,48)[0])
                assert struct.unpack_from('<Q',b,40)[0]==1 and b[48:80]==proof[64:96]
                key=slot,bytes(p[48:128]),bytes(p[368:704])
                branches=int(lo>=0)+int(hi>=0);old_naive_branches+=branches*256
                if key not in groups:
                    groups[key]=len(group_rows)
                    group_rows.append(dict(group=groups[key],slot_id=slot,source=source,target=target,kind='regular_saved_axis1_domain',
                        first_native_panel=index,dependent=dep,selected_branches=branches,maximum_distinct_real_graph_values=16*branches,
                        prepared_packet=str(path),prepared_record=actual,frame_packet='',frame_record=-1))
                row.update(geometry_group=groups[key],prepared_boundary_packet=str(path),prepared_boundary_record=actual)
        else:
            assert r['panel_kind']=='oblique' and len(p)==192
            s,src,tgt,ch,idx,slab,child=struct.unpack_from('<7Q',p);assert (s,src,tgt,idx)==(slot,source,target,index)
            f=frames[slot,ch];frame=record(f['typed_packet'],int(f['record']))
            assert struct.unpack_from('<4Q',frame)==(slot,source,target,ch)
            row.update(oblique_frame_packet=f['typed_packet'],oblique_frame_record=int(f['record']))
            key=slot,'oblique',ch,slab,child,bytes(p[64:128]);classes['oblique']+=1
            if key not in groups:
                groups[key]=len(group_rows)
                # Two is a conservative branch-value cap until actual slab
                # lower/upper identities are bound by the runtime evaluator.
                group_rows.append(dict(group=groups[key],slot_id=slot,source=source,target=target,kind='oblique_same_frame_real_branches',
                    first_native_panel=index,dependent=1,selected_branches=2,maximum_distinct_real_graph_values=32,
                    prepared_packet='',prepared_record=-1,frame_packet=f['typed_packet'],frame_record=int(f['record'])))
            row['geometry_group']=groups[key]
        jobs.append(row)
    assert len(jobs)==84489 and len(slots)==465 and sum(classes.values())==84489
    for pid,p in packets.items():
        if p['path'] in raw:assert sha(Path(p['path']))==p['sha256']
    for name,records in (('panel_jobs.csv',jobs),('real_graph_groups.csv',group_rows)):
        with (dest/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');w.writeheader();w.writerows(records)
    regular_groups=[g for g in group_rows if g['kind']=='regular_saved_axis1_domain']
    regular_values=sum(g['maximum_distinct_real_graph_values'] for g in regular_groups)
    report=dict(scope='exact remaining19target runtime inputs and existing prepared graph-domain links; no geometry/field/physical evaluation, no new arrangement search',
        remaining_targets=sorted(remaining),remaining_slots=465,remaining_panels=84489,
        planned_remaining_native_queries=21629184,maximum_positive_packet_bytes=21629184*1008,
        panel_classes=dict(classes),regular_prepared_graph_groups=len(regular_groups),
        regular_distinct_real_boundary_value_cap=regular_values,
        old_per_point_regular_graph_value_calls=old_naive_branches,
        regular_graph_reuse_reduction_factor=float(old_naive_branches/regular_values),
        oblique_graph_groups=len(group_rows)-len(regular_groups),
        conservative_all_real_boundary_value_cap=sum(g['maximum_distinct_real_graph_values'] for g in group_rows),
        required_method='reuse saved full-real-free-axis implicit domains as initial enclosures; refine only actual Gauss real graph values once per free node, cache across eta panels/nodes',
        preserve_full_event_and_root_uncertainty=True,
        new_canonical_root_or_arrangement_jobs=0,new_field_or_native_queries=0,
        completed_four_native_targets_and_original68526_excluded=True,
        mixed_runtime_geometry_and_error_acquisition_implemented=False,
        full24target_transfer_complete=False,accepted_RK4_evolution=False,long_Kepler_validated=False)
    (dest/'contract.json').write_text(json.dumps(report,indent=2)+'\n')
    (dest/'input.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in dict.fromkeys(inputs)))
    (dest/'source.py').write_bytes(Path(__file__).read_bytes())
    (dest/'source.sha256').write_text(sha(dest/'source.py')+'  '+str(dest/'source.py')+'\n')
    (dest/'output.sha256').write_text(''.join(sha(dest/n)+'  '+str(dest/n)+'\n' for n in ('contract.json','panel_jobs.csv','real_graph_groups.csv')))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()

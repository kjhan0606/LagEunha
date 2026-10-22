"""Adopt actual29 plans/clearance and record evidence for every fixed panel.

No numerical point evaluation. Older plan-only data never becomes a complex
or all12 certificate merely by appearing in this registry.
"""
from collections import Counter,defaultdict
from pathlib import Path
import csv,hashlib,json,struct
from audit_shared_regular_q_transfer_v1 import rows,sha

ROOT=Path(__file__).resolve().parent

def main():
    dest=ROOT/'shared_atlas_adoption_v6';dest.mkdir()
    old=ROOT/'shared_atlas_adoption_v5';plan=ROOT/'shared_pending_atlas_v1_mpi1'
    regular=ROOT/'shared_pending_clearance_v1_mpi1'
    cleared=ROOT/'SHARED_PENDING_REGULAR_CLEARANCE_20261021.json'
    costs=ROOT/'SHARED_PENDING_TYPED_ATLAS_SCOPE_CORRECTION_20261021.json'
    native=json.loads(cleared.read_text());audit=json.loads(costs.read_text())
    review=ROOT/'ASTRA_SHARED_PENDING29_ATLAS_END_REVIEW_20261021.md'
    assert review.is_file() and native['independent_pending_regular_clearance_audit_pass'] and native['all29_new_regular_clearance_complete']
    assert native['new_axis_passes']==101564 and audit['independent_pending29_typed_plan_audit_pass']
    raw_cache={};packets=[];pid={};inputs=[cleared,costs,review,old/'contract.json',old/'effective_slot_plan_references.csv']
    def packet(path,epoch):
        path=Path(path);s=str(path)
        if s not in pid:
            raw=path.read_bytes();h=struct.unpack_from('<8Q',raw)
            assert h[0] in (0x4753415749544e31,0x475341544c415331) and h[1:3]==(1,1) and h[6:]==(24,384) and len(raw)==64+h[3]*h[4]
            pid[s]=len(packets);packets.append(dict(packet_id=pid[s],path=s,record_width=h[3],records=h[4],sha256=hashlib.sha256(raw).hexdigest(),epoch=epoch))
            raw_cache[s]=raw;inputs.append(path)
        return pid[s]
    def record(path,j):
        raw=raw_cache[str(path)];h=struct.unpack_from('<8Q',raw);assert 0<=j<h[4]
        return memoryview(raw)[64+j*h[3]:64+(j+1)*h[3]]
    field={}
    def bind_fields(field_path,axis_path,kind,epoch,width):
        f_id=packet(field_path,epoch);a_id=packet(axis_path,epoch)
        fh=packets[f_id];assert fh['record_width']==width
        for j in range(fh['records']):
            b=record(field_path,j)
            if kind=='regular':slot,source,target,index,axis=struct.unpack_from('<5Q',b)
            else:slot,index,axis=struct.unpack_from('<3Q',b)
            key=slot,kind,index,axis;assert key not in field
            field[key]=(f_id,j,a_id,j)
    bind_fields(regular/'passed_regular_full_fields_rank0.bin',regular/'regular_axis_proofs_rank0.bin','regular','regular_prepared_v1',1680)
    # The historical geometric-v axis file has14 failed records between its
    # passed prefixes; bind its full fields to actual passing axis indices.
    oblique_groups=[(ROOT/'shared_oblique_missing_fields_v1_mpi1','recovered_full_field_points_rank0.bin','recovered_axis_records_rank0.bin','missing800_recovered_v1'),
                    (ROOT/'shared_oblique_v_suffix_v1_mpi1','new_passed_full_field_points_rank0.bin','new_geometric_v_axis_proofs_rank0.bin','geometric_v_prior_v1'),
                    (ROOT/'shared_oblique_collected_suffix_v1_mpi1','new_passed_full_field_points_rank0.bin','new_collected_axis_proofs_rank0.bin','oblique_collected_v1')]
    for folder,fn,an,epoch in oblique_groups:
        fp=folder/fn;ap=folder/an;f_id=packet(fp,epoch);a_id=packet(ap,epoch)
        amap={}
        for j in range(packets[a_id]['records']):
            b=record(ap,j)
            if struct.unpack_from('<Q',b,24)[0]:amap[struct.unpack_from('<3Q',b)]=j
        for j in range(packets[f_id]['records']):
            b=record(fp,j);slot,index,axis=struct.unpack_from('<3Q',b);key=slot,'oblique',index,axis;assert key not in field
            field[key]=(f_id,j,a_id,amap[slot,index,axis])
    assert len(field)==110332
    regular_file=plan/'new_regular_typed_panels_rank0.bin';oblique_file=plan/'unchanged_oblique_panels_rank0.bin'
    reg_id=packet(regular_file,'pending29_typed_v1');obl_id=packet(oblique_file,'adopted_oblique_typed_v1')
    new_by_slot=defaultdict(list)
    for file,kind,p_id in ((regular_file,'regular',reg_id),(oblique_file,'oblique',obl_id)):
        for j in range(packets[p_id]['records']):
            b=record(file,j)
            if kind=='regular':slot,source,index,target=struct.unpack_from('<4Q',b)
            else:slot,source,target,ch,index=struct.unpack_from('<5Q',b)
            new_by_slot[slot].append((kind,index,p_id,j))
    assert len(new_by_slot)==29
    historical_clear={3200,3203,3208,3211}
    old_rows=rows(old/'effective_slot_plan_references.csv');assert len(old_rows)==9216
    extended=['regular_typed_packet','regular_record_first','regular_record_count','oblique_typed_packet','oblique_record_first','oblique_record_count','geometry_density_evidence_status','full_field_evidence_status']
    slot_rows=[];evidence_counts=Counter();pending_axes=[];panel_count=0
    # IDs avoid repeated long path/hash strings in the91168-row evidence table.
    columns=['slot_id','source','target','panel_kind','native_panel','typed_packet_id','typed_record',
        'domain_status','complex_axes_status','full_field_status','axis0_field_packet','axis0_field_record',
        'axis0_proof_packet','axis0_proof_record','axis1_field_packet','axis1_field_record','axis1_proof_packet','axis1_proof_record',
        'thermal_context_status','all12_majorant_status','quadrature_remainder_status','endpoint_error_status',
        'native_forward_sum_final_status','integral_status']
    with (dest/'panel_evidence_registry.csv').open('w',newline='') as pf:
        writer=csv.writer(pf);writer.writerow(columns)
        def emit(r,kind,index,p_id,j,domain,complex_status,field_status):
            nonlocal panel_count
            slot=int(r['slot_id']);source=int(r['source']);target=int(r['target']);links=[]
            for axis in (0,1):
                k=slot,kind,index,axis
                if field_status=='stored_pair':assert k in field;links.extend(field[k])
                else:links.extend((-1,-1,-1,-1))
            published=kind=='retained_group'
            writer.writerow((slot,source,target,kind,index,p_id,j,domain,complex_status,field_status,*links,
                'retained_group' if published else 'binding_pending',
                *(['retained_group']*4 if published else ['pending']*4),
                'retained_original_target_packet' if published else 'not_evaluated'))
            evidence_counts[field_status]+=1;panel_count+=1
            if field_status in ('missing_snapshot','unbound'):
                for axis in (0,1):pending_axes.append((slot,source,target,kind,index,axis,p_id,j,complex_status))
        for before in old_rows:
            r=dict(before);slot=int(r['slot_id']);n=int(r['active_panels'])
            for name in extended:r[name]=''
            if slot in new_by_slot:
                r.update(atlas_complete='1',active_panels=str(len(new_by_slot[slot])),complete_planned_queries=str(256*len(new_by_slot[slot])),provenance='adopted_pending29_same_function_real_plan_and_geometry_density_clearance',typed_packet='HYBRID_SEE_EXPLICIT_PACKET_REFERENCES',typed_record_first='-1',typed_record_count=str(len(new_by_slot[slot])),analytic_ready='0',required_reference_radius_evaluator='origin_quadratic_or_original_same_function',required_angular_proof_evaluator='prepared_weighted_regular_or_oblique_collected',radius_branch_proof_reference=str(regular/'passed_regular_full_fields_rank0.bin'))
                for kind,index,p_id,j in new_by_slot[slot]:emit(r,kind,index,p_id,j,'adopted_actual_chart_coverage','passed_both_axes','stored_pair')
                for kind,filename in (('regular',regular_file),('oblique',oblique_file)):
                    v=[t for t in new_by_slot[slot] if t[0]==kind]
                    if v:r[kind+'_typed_packet']=str(filename);r[kind+'_record_first']=str(v[0][3]);r[kind+'_record_count']=str(len(v));assert [x[3] for x in v]==list(range(v[0][3],v[0][3]+len(v)))
                r['geometry_density_evidence_status']='passed_both_axes';r['full_field_evidence_status']='stored_pair'
            elif n:
                if r['target']=='68526':
                    for j in range(n):emit(r,'retained_group',j,-1,-1,'retained_target_publication','retained_group_not_individual_axis','retained_group')
                    r['geometry_density_evidence_status']='retained_original_target_publication';r['full_field_evidence_status']='retained_group'
                else:
                    path=Path(r['typed_packet']);p_id=packet(path,'inherited_v5_plan');first=int(r['typed_record_first']);assert int(r['typed_record_count'])==n
                    for j in range(first,first+n):
                        b=record(path,j);actual_slot,source,index,target=struct.unpack_from('<4Q',b);assert (actual_slot,source,target)==(slot,int(r['source']),int(r['target']))
                        emit(r,'regular',index,p_id,j,'inherited_actual_plan',
                            'historical_pass_no_full_field' if slot in historical_clear else 'not_bound_by_plan_registry',
                            'missing_snapshot' if slot in historical_clear else 'unbound')
                    r['geometry_density_evidence_status']='historical_source133_pass' if slot in historical_clear else 'not_decided_by_plan'
                    r['full_field_evidence_status']='missing_snapshot' if slot in historical_clear else 'unbound'
            else:r['geometry_density_evidence_status']='retained_proved_empty';r['full_field_evidence_status']='not_required_empty'
            assert r['atlas_complete']=='1';slot_rows.append(r)
    assert panel_count==91168 and evidence_counts=={'stored_pair':55166,'retained_group':1047,'missing_snapshot':3976,'unbound':30979}
    assert len(pending_axes)==69910
    with (dest/'effective_slot_plan_references.csv').open('w',newline='') as f:
        out=csv.DictWriter(f,fieldnames=list(old_rows[0])+extended);out.writeheader();out.writerows(slot_rows)
    with (dest/'unbound_complex_axis_inputs.csv').open('w',newline='') as f:
        out=csv.writer(f);out.writerow(('slot_id','source','target','panel_kind','native_panel','axis','typed_packet_id','typed_record','prior_clearance_status'));out.writerows(pending_axes)
    # This is a binding/recovery candidate list, NOT authority to blindly
    # replay69910 queries if an existing stronger certificate can be linked.
    packet_fields=list(packets[0])
    with (dest/'typed_packet_registry.csv').open('w',newline='') as f:
        out=csv.DictWriter(f,fieldnames=packet_fields);out.writeheader();out.writerows(packets)
    target_cost=[]
    for r in audit['target_cost']:
        row=dict(r);target=int(row['target']);n=sum(int(s['active_panels']) for s in slot_rows if int(s['target'])==target)
        assert n==row['complete_candidate_plan_panels'] and n*256==row['planned_positive_queries']
        row['geometry_density_clearance_status']='retained_published_original_target' if target==68526 else 'panel_registry_required_not_inferred_from_cost'
        row['all12_analytic_remainders_complete']=target==68526
        target_cost.append(row)
    # Refer to the retained whole-target publication once. Never invent
    # per-source or per-panel integrals by dividing a group result.
    retained_target=ROOT/'shared_missing_source_v1_mpi1/caller_target_packet.bin'
    assert retained_target.is_file();inputs.append(retained_target)
    result=dict(scope='all9216 source-target typed plans complete;29 formerly missing atlases newly adopted with their geometry/density evidence; per-panel evidence remains explicit, no all24 complex/all12 inference',
        previous_registry=str(old/'contract.json'),completed_END_review=str(review),review_sha256=sha(review),
        conditional_clearance_audit_satisfied_by=str(cleared),all_target_actual_atlas_complete=True,
        active_new_panels=90121,retained_original_target_panels=1047,total_panels=91168,
        remaining_incomplete_atlas_slots=0,all24_candidate_costs_known=True,
        maximum_target_queries=4496896,authorized_target_query_limit=13074432,
        newly_adopted_slots=sorted(new_by_slot),new_geometry_density_regular_axis_passes=101564,
        reused_oblique_full_fields=8768,existing800_oblique_recovery_already_complete=True,
        stored_new29_full_fields=110332,panel_full_field_status_counts=dict(evidence_counts),
        unbound_or_missing_complex_axis_records=69910,unbound_list_is_not_a_blind_replay_instruction=True,
        source133_historical_pass_axes_needing_missing_field_binding=7952,
        baseline_plan_axes_without_bound_complex_evidence=61958,
        retained_original_target_publication=dict(target=68526,path=str(retained_target),sha256=sha(retained_target),scope='whole target only; grouped source data not split into individual integrals'),
        all_target_geometry_density_clearance_complete=False,all_target_analytic_certificates_complete=False,
        full_quadrature_allowed=False,new_physical_quadrature_queries=0,accepted_RK4_evolution=False,long_Kepler_validated=False,target_cost=target_cost)
    (dest/'contract.json').write_text(json.dumps(result,indent=2)+'\n')
    inputs.extend((ROOT/'SHARED_PENDING29_CHART_ASSEMBLY_20261021.json',ROOT/'SHARED_OBLIQUE_MISSING_FIELDS_20261021.json',ROOT/'shared_oblique_all12_inputs_v1/contract.json'))
    (dest/'input.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in dict.fromkeys(inputs)))
    (dest/'source.sha256').write_text(sha(Path(__file__))+'  '+str(Path(__file__).resolve())+'\n')
    names=('contract.json','effective_slot_plan_references.csv','panel_evidence_registry.csv','unbound_complex_axis_inputs.csv','typed_packet_registry.csv')
    (dest/'output.sha256').write_text(''.join(sha(dest/n)+'  '+str(dest/n)+'\n' for n in names))
    print(json.dumps({k:v for k,v in result.items() if k!='target_cost'},indent=2))

if __name__=='__main__':main()

"""Bind audited baseline/suffix fields to v6 panels without replay or replanning."""
from collections import Counter
from pathlib import Path
import csv,json,struct
from audit_shared_regular_q_transfer_v1 import rows,sha

ROOT=Path(__file__).resolve().parent

def main():
    old=ROOT/'shared_atlas_adoption_v6';dest=ROOT/'shared_atlas_adoption_v7';dest.mkdir()
    baseline=ROOT/'shared_baseline_fields_audit_v1/contract.json'
    suffix=ROOT/'shared_baseline_density_suffix_audit_v1/contract.json'
    density=ROOT/'shared_density_witness_audit_v1/contract.json'
    b=json.loads(baseline.read_text());s=json.loads(suffix.read_text());d=json.loads(density.read_text())
    assert b['independent_baseline_field_audit_pass'] and b['new_axis_passes']==58646
    assert s['independent_baseline_density_suffix_audit_pass'] and s['new_axis_passes']==11264
    assert d['independent_same_radius_density_accounting_pass']
    inputs=[baseline,suffix,density]+[old/n for n in ('contract.json','panel_evidence_registry.csv','typed_packet_registry.csv','effective_slot_plan_references.csv')]
    packets=rows(old/'typed_packet_registry.csv');field={}
    for name,epoch,required in (('shared_baseline_fields_v2_mpi1','baseline_regular_prepared_v2',58646),
                                 ('shared_baseline_density_suffix_v2_mpi1','baseline_same_density_two_sided_v2',11264)):
        run=ROOT/name;fp=run/'passed_regular_full_fields_rank0.bin';ap=run/'regular_axis_proofs_rank0.bin'
        native=json.loads((run/'contract.json').read_text());assert native['saved_full_fields']==required
        raw=[];ids=[]
        for path,width in ((fp,1680),(ap,624)):
            data=path.read_bytes();h=struct.unpack_from('<8Q',data)
            assert h[:3]==(0x4753415749544e31,1,1) and h[3]==width and len(data)==64+h[3]*h[4]
            raw.append(data);ids.append(len(packets))
            packets.append(dict(packet_id=len(packets),path=str(path),record_width=width,records=h[4],sha256=sha(path),epoch=epoch));inputs.append(path)
        proofs={}
        for j in range((len(raw[1])-64)//624):
            at=64+j*624;key=struct.unpack_from('<5Q',raw[1],at)
            if struct.unpack_from('<Q',raw[1],at+40)[0]:assert key not in proofs;proofs[key]=j
        assert len(proofs)==required
        for j in range(required):
            key=struct.unpack_from('<5Q',raw[0],64+j*1680);slot,source,target,index,axis=key
            assert key in proofs and (slot,index,axis) not in field
            field[slot,index,axis]=(ids[0],j,ids[1],proofs[key],source,target,epoch)
    assert len(field)==69910
    counts=Counter();bound=0
    with (old/'panel_evidence_registry.csv').open(newline='') as f,(dest/'panel_evidence_registry.csv').open('w',newline='') as out:
        reader=csv.DictReader(f);writer=csv.DictWriter(out,fieldnames=reader.fieldnames);writer.writeheader()
        for row in reader:
            if row['full_field_status'] in ('unbound','missing_snapshot'):
                slot,index=int(row['slot_id']),int(row['native_panel'])
                for axis in (0,1):
                    x=field[slot,index,axis];assert (x[4],x[5])==(int(row['source']),int(row['target']))
                    for key,value in zip(('field_packet','field_record','proof_packet','proof_record'),x[:4]):row[f'axis{axis}_{key}']=str(value)
                    bound+=1
                row['complex_axes_status']='passed_both_axes'
                row['full_field_status']='stored_pair'
            counts[row['full_field_status']]+=1;writer.writerow(row)
    assert bound==69910 and counts=={'stored_pair':90121,'retained_group':1047}
    slot_rows=rows(old/'effective_slot_plan_references.csv')
    for r in slot_rows:
        if r['full_field_evidence_status'] in ('unbound','missing_snapshot'):
            r['full_field_evidence_status']='stored_pair';r['geometry_density_evidence_status']='passed_both_axes'
    assert not any(r['full_field_evidence_status'] in ('unbound','missing_snapshot') for r in slot_rows)
    with (dest/'effective_slot_plan_references.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(slot_rows[0]));w.writeheader();w.writerows(slot_rows)
    with (dest/'typed_packet_registry.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(packets[0]));w.writeheader();w.writerows(packets)
    c=json.loads((old/'contract.json').read_text())
    c.update(scope='same complete v6 actual integration plans, now all active panels have audited two-axis geometry/density evidence and full sd::Point/radius/measure data; not all12, quadrature or evolution',
       previous_registry=str(old/'contract.json'),baseline_clearance_audit=str(baseline),
       same_density_suffix_audit=str(suffix),same_density_exponential_identity_audit=str(density),
       panel_full_field_status_counts=dict(counts),stored_active_full_fields=180242,
       newly_bound_baseline_full_fields=69910,baseline_new_axis_checks=58710,
       density_suffix_new_axis_checks=11264,first_density_failure_disks_rechecked_without_field_queries=64,
       required_reference_density_evaluator=str(ROOT/'shared_baseline_density_suffix_v2_mpi1/material_shared_logistic_density.h'),
       required_reference_density_evaluator_sha256=sha(ROOT/'shared_baseline_density_suffix_v2_mpi1/material_shared_logistic_density.h'),
       remaining_incomplete_atlas_slots=0,unbound_or_missing_complex_axis_records=0,
       source133_historical_pass_axes_needing_missing_field_binding=0,
       baseline_plan_axes_without_bound_complex_evidence=0,
       all_target_geometry_density_clearance_complete=True,
       all_target_analytic_certificates_complete=False,full_quadrature_allowed=False,
       frozen_thermal_context_still_requires_explicit_binding=True,
       general_evolved_heat_restart_supported=False,new_physical_quadrature_queries=0,
       accepted_RK4_evolution=False,long_Kepler_validated=False)
    for row in c['target_cost']:
        if row['target']!=68526:row['geometry_density_clearance_status']='all_actual_panel_axes_audited_and_full_fields_stored'
    (dest/'contract.json').write_text(json.dumps(c,indent=2)+'\n')
    (dest/'input.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in dict.fromkeys(inputs)))
    (dest/'source.py').write_bytes(Path(__file__).read_bytes())
    (dest/'source.sha256').write_text(sha(dest/'source.py')+'  '+str(dest/'source.py')+'\n')
    names=('contract.json','effective_slot_plan_references.csv','panel_evidence_registry.csv','typed_packet_registry.csv')
    (dest/'output.sha256').write_text(''.join(sha(dest/n)+'  '+str(dest/n)+'\n' for n in names))
    print(json.dumps(dict(bound_new_full_fields=bound,total_active_full_fields=180242,
       panel_status_counts=dict(counts),remaining_missing_axis_inputs=0,
       all_target_geometry_density_clearance_complete=True,
       all_target_analytic_certificates_complete=False,full_quadrature_allowed=False)))

if __name__=='__main__':main()

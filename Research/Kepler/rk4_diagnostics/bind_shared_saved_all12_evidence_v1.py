"""Bind the audited all12 records back to exact saved measures and J values."""
from pathlib import Path
import json,struct
from audit_shared_regular_q_transfer_v1 import rows,sha

ROOT=Path(__file__).resolve().parent

def main():
    dest=ROOT/'shared_saved_all12_evidence_v1';dest.mkdir()
    registry=ROOT/'shared_atlas_adoption_v7';run=ROOT/'shared_saved_all12_v1_mpi1'
    audit=ROOT/'shared_saved_all12_audit_v1/contract.json'
    a=json.loads(audit.read_text());assert a['independent_saved_all12_analytic_accounting_pass']
    thermal=ROOT/'shared_all12_control_context_v1/contract.json'
    t=json.loads(thermal.read_text());assert t['frozen_CONTROL_thermal_context_bound']
    context={int(r['source']):r for r in rows(ROOT/'shared_all12_control_context_v1/source_thermal_context_bindings.csv')}
    packets={int(p['packet_id']):p for p in rows(registry/'typed_packet_registry.csv')}
    cache={};inputs=[audit,thermal,registry/'contract.json',registry/'typed_packet_registry.csv',registry/'panel_evidence_registry.csv',run/'saved_field_all12_panel_proofs_rank0.bin']
    def packet(pid):
        if pid not in cache:
            p=packets[pid];path=Path(p['path']);data=path.read_bytes();h=struct.unpack_from('<8Q',data)
            assert h[3]==int(p['record_width']) and h[4]==int(p['records']) and len(data)==64+h[3]*h[4]
            assert sha(path)==p['sha256'];cache[pid]=(data,h[3],h[4]);inputs.append(path)
        return cache[pid]
    def record(pid,index):
        data,width,count=packet(pid);assert 0<=index<count
        return memoryview(data)[64+width*index:64+width*(index+1)]
    proof=(run/'saved_field_all12_panel_proofs_rank0.bin').read_bytes();index=0;measure_links=0
    for row in rows(registry/'panel_evidence_registry.csv'):
        if row['panel_kind']=='retained_group':continue
        p=memoryview(proof)[64+1664*index:64+1664*(index+1)];index+=1
        source=int(row['source']);assert source in context
        typed=record(int(row['typed_packet_id']),int(row['typed_record']))
        assert struct.unpack_from('<Q',p,32)[0]==(typed[192] if row['panel_kind']=='regular' else 0)
        for axis in (0,1):
            f=record(int(row[f'axis{axis}_field_packet']),int(row[f'axis{axis}_field_record']))
            if row['panel_kind']=='regular':J,measure=432,1472
            else:J,measure=416,1456
            assert p[832+32*axis:864+32*axis]==f[measure:measure+32]
            volume=64+384*axis+32*6
            assert p[volume:volume+32]==f[J:J+32]
            measure_links+=1
    assert index==90121 and measure_links==180242
    summary=dict(scope='audited all12 proof records connected to original plan whole flags, stored reference measures/current J and384 immutable CONTROL/source EOS contexts; preserves native analytic/libm premises',
       adopted_registry=str(registry/'contract.json'),all12_analytic_accounting_audit=str(audit),
       panel_plan_subtraction_method_links=90121,axis_reference_measure_and_current_J_links=180242,
       source_thermal_contexts=384,all_new23_target_all12_analytic_subtotals_within_original_tau=True,
       original68526_publication_retained=True,
       maximum_analytic_subtotal_fraction_of_original_tau=a['worst_analytic_subtotal_fractions'][0],
       all12_panel_remainder_inputs_complete=True,
       endpoint_domain_error_complete=False,native_forward_sum_final_complete=False,
       total_target_error_budget_complete=False,full_quadrature_allowed=False,
       general_evolved_heat_restart_supported=False,heat_reset_permitted=False,
       new_complex_field_queries=0,new_physical_quadrature_queries=0,
       accepted_RK4_evolution=False,long_Kepler_validated=False,
       completed_END_review=str(ROOT/'ASTRA_SHARED_BASELINE_DENSITY_END_REVIEW_20261021.md'))
    (dest/'contract.json').write_text(json.dumps(summary,indent=2)+'\n')
    (dest/'input.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in dict.fromkeys(inputs)))
    (dest/'source.py').write_bytes(Path(__file__).read_bytes())
    (dest/'source.sha256').write_text(sha(dest/'source.py')+'  '+str(dest/'source.py')+'\n')
    (dest/'output.sha256').write_text(sha(dest/'contract.json')+'  '+str(dest/'contract.json')+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()

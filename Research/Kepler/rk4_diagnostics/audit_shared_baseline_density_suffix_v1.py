"""Only the64failed slots: saved exact panels and full-field accounting."""
from pathlib import Path
import contextlib,io,json,struct
from audit_shared_atlas_witness_v1 import wire
from audit_shared_pending29_projection_v1 import verify
from audit_shared_regular_q_transfer_v1 import rows
import audit_shared_pending_clearance_v1 as regular

ROOT=Path(__file__).resolve().parent

def main():
    name='shared_baseline_density_suffix_v2_mpi1';run=ROOT/name;c=verify(run)
    old=ROOT/'shared_baseline_fields_v2_mpi1'
    _,axes=wire(old/'regular_axis_proofs_rank0.bin')
    failed={struct.unpack_from('<Q',a)[0] for a in axes if not struct.unpack_from('<Q',a,40)[0]}
    assert len(failed)==64
    assert not any(struct.unpack_from('<Q',a,40)[0] for a in axes if struct.unpack_from('<Q',a)[0] in failed)
    _,old_panels=wire(old/'unchanged_baseline_panels_rank0.bin')
    selected=[p for p in old_panels if struct.unpack_from('<Q',p)[0] in failed]
    _,panels=wire(run/'unchanged_baseline_panels_rank0.bin')
    assert panels==selected and len(panels)==5632
    density=json.loads((ROOT/'shared_density_witness_audit_v1/contract.json').read_text())
    assert density['independent_same_radius_density_accounting_pass'] and density['same_density_disk_passes']==64
    output=io.StringIO()
    with contextlib.redirect_stdout(output):
        regular.main(name,name,'unchanged_baseline_panels_rank0.bin',5632,64)
    report=json.loads(output.getvalue())
    _,proof=wire(run/'selected_derivative_proofs_rank0.bin')
    assert len(proof)==c['selected_derivative_records']==0
    assert c['new_complex_axis_evaluations']==c['new_complex_axis_passes']==c['saved_full_fields']==11264
    assert c['completed_regular_slot_clearances']==64 and not c['density_profile_changed']
    report.update(independent_baseline_density_suffix_audit_pass=True,
                  unchanged_original_panel_bytes=5632,previous_successful_field_replays=0,
                  all64_previously_failed_slots_complete=True,
                  same_density_two_sided_logistic_proof=str(ROOT/'shared_density_witness_audit_v1/contract.json'),
                  primitive_logistic_domain_checked_for_every_new_field=True,
                  scope='exact unproved64slot plan/full-field/cache/branch accounting, same two-sided logistic primitive for every axis; native analytic/libm premises retained, not all12 or trajectory')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()

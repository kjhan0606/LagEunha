"""Exact missing800 selection and complete14-subchart field packet binding."""
from pathlib import Path
import json
import struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball,add,mul,neg,encloses,verify
from audit_shared_regular_q_transfer_v1 import sha
from audit_shared_oblique_v_suffix_v1 import audit_domains

ROOT=Path(__file__).resolve().parent

def main():
    run=ROOT/'shared_oblique_missing_fields_v1_mpi1'
    prior=ROOT/'shared_oblique_v_suffix_v1_mpi1'
    adopted=ROOT/'shared_oblique_collected_suffix_v1_mpi1'
    c=verify(run)
    def read(name,width):
        h,records=wire(run/name)
        assert h[3]==width
        return records
    selected=read('selected_original800_inputs_rank0.bin',656)
    assert selected==wire(prior/'retained800_axis_proofs_rank0.bin')[1]
    axes=read('recovered_axis_records_rank0.bin',656)
    points=read('recovered_full_field_points_rank0.bin',1488)
    domains=read('recovery_geometric_domains_rank0.bin',880)
    _,frames=wire(adopted/'unchanged_frames_rank0.bin')
    _,panels=wire(adopted/'unchanged_candidate_panels_rank0.bin')
    _,slabs=wire(adopted/'unchanged_actual_slabs_rank0.bin')
    fm={struct.unpack_from('<Q',b)[0]:b for b in frames}
    pm={(struct.unpack_from('<Q',b)[0],struct.unpack_from('<Q',b,32)[0]):b for b in panels}
    sm={(struct.unpack_from('<Q',b)[0],struct.unpack_from('<Q',b,32)[0]):b for b in slabs}
    key=lambda b:struct.unpack_from('<3Q',b)
    assert len(selected)==len(axes)==len(points)==800
    assert {key(b) for b in axes}=={key(b) for b in selected}
    am={key(b):b for b in axes}
    for old,b,p in zip(selected,axes,points):
        slot,index,axis=key(b)
        assert key(b)==key(old)==key(p) and old[32:96]==b[32:96]
        assert struct.unpack_from('<Q',b,24)[0]==1
        assert b[256:320]==p[1424:1488] and b[224:256]==p[416:448]
        assert ball(b,224)[0]>ball(b,224)[1]>=0
        assert 0<=quad(b,352)<1 and 0<=quad(b,368)<1
        assert 0<=quad(b,400)<=quad(b,416)
        panel=pm[slot,index]
        measure=mul(mul(mul(ball(fm[slot],576),ball(b,192)),
                        add(ball(panel,96),neg(ball(panel,64)))),
                    add(ball(panel,160),neg(ball(panel,128))))
        encloses(ball(b,288),measure)
    state_path=ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin'
    assert sha(state_path)=='74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005'
    audit_domains(domains,fm,state_path.read_bytes())
    for z in domains:
        slot,source,target,branch,index,axis=struct.unpack_from('<6Q',z)
        assert z[64:96]==am[slot,index,axis][32:64]
        slab=sm[slot,struct.unpack_from('<Q',pm[slot,index],40)[0]]
        assert branch in struct.unpack_from('<2q',slab,64)
    linked={key(p):(str(run/'recovered_full_field_points_rank0.bin'),i)
            for i,p in enumerate(points)}
    for owner in (prior,adopted):
        file=owner/'new_passed_full_field_points_rank0.bin'
        _,stored=wire(file)
        for i,p in enumerate(stored):
            assert key(p) not in linked
            linked[key(p)]=str(file),i
    assert len(linked)==8768
    assert set(linked)=={(slot,index,axis) for slot,index in pm for axis in (0,1)}
    assert c['complete800_recovery'] and c['new_required_complex_field_queries']==800
    assert c['unselected_field_queries']==c['existing7968_full_fields_recomputed']==c['new_physical_quadrature_queries']==c['new_inverse_solves']==c['new_canonical_root_jobs']==0
    assert c['state_unchanged'] and not c['full_quadrature_allowed']
    print(json.dumps(dict(scope='exact800 originally missing complex inputs/field/J/radius/measure and new geometric-domain bounds; disjoint800+2005+5963 complete bindings; retained native analytic premises, not all12 or time evolution',
        run=str(run),independent_missing800_field_audit_pass=True,
        recovered_required_fields=800,retained_full_fields=7968,total_full_field_points=8768,
        complete14_oblique_full_field_bindings=True,new_boundary_neighborhood_records=len(domains),
        field_packets=[dict(path=str(file),sha256=sha(file)) for file in (
            run/'recovered_full_field_points_rank0.bin',prior/'new_passed_full_field_points_rank0.bin',adopted/'new_passed_full_field_points_rank0.bin')],
        remaining_incomplete_atlas_slots=29,all12_analytic_remainders_complete=False,
        full_quadrature_allowed=False,long_time_stability_demonstrated=False),indent=2))

if __name__=='__main__':main()

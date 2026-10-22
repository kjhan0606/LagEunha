"""Complete unchanged14-chart complex clearance accounting, not RK4 physics."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import json
import struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball,add,mul,neg,encloses,verify
from audit_shared_regular_q_transfer_v1 import rows,sha
from audit_shared_oblique_v_suffix_v1 import audit_domains

ROOT=Path(__file__).resolve().parent

def main():
    run=ROOT/'shared_oblique_collected_suffix_v1_mpi1'
    prior=ROOT/'shared_oblique_v_suffix_v1_mpi1'
    original=ROOT/'shared_oblique_panels_v1_mpi1'
    c=verify(run)
    def read(name,width):
        h,records=wire(run/name)
        assert h[3]==width
        return records
    frames=read('unchanged_frames_rank0.bin',2064)
    panels=read('unchanged_candidate_panels_rank0.bin',192)
    slabs=read('unchanged_actual_slabs_rank0.bin',304)
    assert frames==wire(prior/'unchanged_frames_rank0.bin')[1]
    assert panels==wire(prior/'unchanged_candidate_panels_rank0.bin')[1]
    assert slabs==wire(prior/'unchanged_actual_slabs_rank0.bin')[1]
    linear=read('same_chart_linear_coordinates_rank0.bin',144)
    assert linear==wire(original/'same_chart_linear_coordinates_rank0.bin')[1]
    retained=read('retained2805_axis_proofs_rank0.bin',656)
    axes=read('new_collected_axis_proofs_rank0.bin',656)
    points=read('new_passed_full_field_points_rank0.bin',1488)
    domains=read('geometric_v_prior_domains_rank0.bin',880)
    old=wire(prior/'retained800_axis_proofs_rank0.bin')[1]
    old += [b for b in wire(prior/'new_geometric_v_axis_proofs_rank0.bin')[1]
            if struct.unpack_from('<Q',b,24)[0]]
    def key(b): return struct.unpack_from('<3Q',b)
    assert len(retained)==len(old)==2805 and {key(b):b for b in retained}=={key(b):b for b in old}
    assert len(axes)==len(points)==c['new_complex_axis_passes']==5963
    assert len(domains)==c['geometric_prior_domain_records']==348
    assert all(struct.unpack_from('<Q',b,24)[0] for b in axes)
    all_axis={key(b):b for b in retained+axes}
    assert len(all_axis)==2*len(panels)==8768
    fm={struct.unpack_from('<Q',b)[0]:b for b in frames}
    pb=defaultdict(list)
    for b in panels: pb[struct.unpack_from('<Q',b)[0]].append(b)
    kept={key(b) for b in retained}
    ab=defaultdict(list)
    for b,p in zip(axes,points):
        slot,index,axis=key(b)
        assert key(b) not in kept and b[:24]==p[:24] and b[256:320]==p[1424:1488]
        ab[slot].append(b)
        panel=pb[slot][index]
        U=(F(1,2),F(1,2) if axis else F(5,4))
        E=(F(1,2),F(5,4) if axis else F(1,2))
        encloses(ball(b,32),add(ball(panel,64),mul(add(ball(panel,96),neg(ball(panel,64))),U)))
        encloses(ball(b,64),add(ball(panel,128),mul(add(ball(panel,160),neg(ball(panel,128))),E)))
        assert ball(b,224)[0]>ball(b,224)[1]>=0
        assert 0<=quad(b,352)<1 and 0<=quad(b,368)<1
        assert 0<=quad(b,400)<=quad(b,416)
        measure=mul(mul(mul(ball(fm[slot],576),ball(b,192)),
                        add(ball(panel,96),neg(ball(panel,64)))),
                    add(ball(panel,160),neg(ball(panel,128))))
        encloses(ball(b,288),measure)
    ledger={int(r['slot_id']):r for r in rows(run/'collected_suffix_clearance.csv')}
    assert set(ledger)==set(pb)==set(fm) and len(fm)==14
    for slot,pan in pb.items():
        expected=[(slot,j,a) for j in range(len(pan)) for a in (0,1)]
        assert all(k in all_axis for k in expected)
        assert [key(b) for b in ab[slot]]==[k for k in expected if k not in kept]
        r=ledger[slot]
        assert r['clearance_complete']=='1' and int(r['new_axis_checks'])==len(ab[slot])
        assert int(r['retained_passing_axes'])+len(ab[slot])==2*len(pan)
    state_path=ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin'
    assert sha(state_path)=='74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005'
    state=state_path.read_bytes()
    audit_domains(domains,fm,state)
    sm={(struct.unpack_from('<Q',b)[0],struct.unpack_from('<Q',b,32)[0]):b for b in slabs}
    for z in domains:
        slot,source,target,branch,index,axis=struct.unpack_from('<6Q',z)
        assert z[64:96]==all_axis[slot,index,axis][32:64]
        panel=pb[slot][index]
        slab=sm[slot,struct.unpack_from('<Q',panel,40)[0]]
        assert branch in struct.unpack_from('<2q',slab,64)
    assert c['passing_axis_replays']==c['new_inverse_solves']==c['new_canonical_root_jobs']==c['new_physical_quadrature_queries']==0
    assert c['completed_oblique_chart_geometry_density_clearance']==14 and c['state_unchanged']
    assert not c['full_quadrature_allowed'] and c['remaining_incomplete_atlas_slots']==29
    # ABI changes are confined to a private helper; none of the read historical
    # typed packets embeds AffineAngular. Historical ray-trace ABI is excluded.
    print(json.dumps(dict(scope='same4384panels/all8768axes: exact prefix preservation, complete suffix ownership, field snapshot links, J/measure and348 geometric pole/binomial/Lipschitz domains; native field enclosure premises retained; not physical evolution',
        run=str(run),independent_collected_oblique_suffix_audit_pass=True,
        retained_passing_axes=2805,new_passing_axes=5963,total_passing_axes=8768,
        unchanged_panels=4384,completed_oblique_chart_geometry_density_clearance=14,
        complete_geometric_prior_domains=348,saved_full_field_points=2005+5963,
        original800_full_field_values_still_needed_for_moment_bounds=True,
        remaining_incomplete_atlas_slots=29,all12_analytic_remainders_complete=False,
        full_quadrature_allowed=False,long_time_stability_demonstrated=False),indent=2))

if __name__=='__main__': main()
